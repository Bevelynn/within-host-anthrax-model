import numpy as np
from scipy.integrate import odeint
import pandas

#timestep for tau leaping
tau_dt = 0.01

##### load dose distribution data #####
dose_distrbution = pandas.read_csv('Sverdlovsk_dose_distribution.csv',header=0)
log_doses = dose_distrbution['log_10(dose)'].to_numpy()
cum_probs = dose_distrbution['cumulative probability'].to_numpy()

#function used to sample a dose from the distribution
def sample_dose():
    urv = np.random.uniform()
    #the following defines the inverse of the cumulative distribution function for the dose distribution
    #returns i such that cum_probs[i-1]<= urv < cum_probs[i]
    idx_2 = np.searchsorted(cum_probs, urv, side = 'right')
    idx_1 = idx_2-1
    log_dose = log_doses[idx_1] + (log_doses[idx_2]-log_doses[idx_1])*(urv - cum_probs[idx_1])/(cum_probs[idx_2]-cum_probs[idx_1])
    return pow(10,log_dose)

#####################################
'''functions to simulate the model'''
#####################################

# deterministic part of the model to find the time (in hours) between
# reaching the threshold bacterial level at which the model switches
# to deterministic and the time at which the expected amount of bacteria reaches 10^10
def derivative(x, t, params):
    #variables are:
    #spores in airways that will lead to infected cells that rupture
    #infected cells in stage 1-3
    #bacterial CFU
    
    S, P1, P2, P3, B = x
    
    dS = -params['rho']*S
    dP1 = params['rho']*S - params['delta']*P1
    dP2 = params['delta']*(P1-P2)
    dP3 = params['delta']*(P2-P3)
    dB = params['delta']*params['R']*P3 + (params['lamb']-params['mu'])*B
    
    return np.hstack((dS, dP1, dP2, dP3, dB))

#solve the model differential eqns
def get_ODE_solution(initial_cond, det_times, params):
    
    sol, infodict = odeint(derivative, initial_cond, det_times, args=(params,), full_output = 1)
    
    B = sol[:,4]

    return B


#functions for the different events that can happen in the stochastic process
def M_eat_S_n(n, S, P1, P2, P3, B, time, params):
    """ macrophages phagocytose n spores """
    n = min(n,S)
    return S-n, P1+n, P2, P3, B

def M_1_2_n(n, S, P1, P2, P3, B, time, params):
    """ n macrophages go from stage 1 to stage 2 """
    n = min(n,P1)
    return S, P1-n, P2+n, P3, B

def M_2_3_n(n, S, P1, P2, P3, B, time, params):
    """ n macrophages go from stage 2 to stage 3 """
    n = min(n,P2)
    return S, P1, P2-n, P3+n, B

def M_rupt_n(n, S, P1, P2, P3, B, time, params):
   """ n host macrophages infected with spores rupture and release bacteria"""
   for _ in range(min(n,P3)):
       ruptsize=np.random.geometric(1/params['R'])
       P3-=1
       B+=ruptsize
   return S, P1, P2, P3, B

def B_birth_n(n, S, P1, P2, P3, B, time, params):
    """ n extracellular bacteria replicate """
    return S, P1, P2, P3, B+n

def B_death_n(n, S, P1, P2, P3, B, time, params):
    """ n extracellular bacteria die """
    return S, P1, P2, P3, max(0,B-n)

events={0:M_eat_S_n, 1:M_1_2_n, 2:M_2_3_n, 3:M_rupt_n, 4:B_birth_n, 5:B_death_n}

#rate of each event
def makerates(S, P1, P2, P3, B, params):
   myrates = []
   myrates.append(params['rho']*S)                            # phagocytosis of spores
   myrates.append(params['delta']*P1)                         # next stage for infected cells
   myrates.append(params['delta']*P2)                         # next stage for infected cells
   myrates.append(params['delta']*P3)                         # next stage for infected cells
   myrates.append(params['lamb']*B)                           # replication of extrac bacteria
   myrates.append(params['mu']*B)                             # death of extrac bacteria
   return myrates

#one step of the Gillespie or tau-leaping algorithm
def onestep(t, S, P1, P2, P3, B, params):
    rates = makerates(S, P1, P2, P3, B, params)
    totalrate = sum(rates)
    #if the average time to do one event is very small (e.g less than 0.001) then do tau leaping
    if 1/totalrate<0.001:
        useGillespie = False
    else:
        useGillespie = True
    if useGillespie:
       urv = np.random.uniform()
       totalrate = sum(rates)
       i = np.searchsorted(np.cumsum(rates),urv*totalrate)
       S, P1, P2, P3, B = events.get(i)(1,S, P1, P2, P3, B, t, params)
       dt = -np.log(np.random.uniform())/totalrate 
    else:
       dt = tau_dt
       mevs = [dt*rate for rate in rates]  # mean numbers of events
       #choose number of events from poisson dist with lambda=mean number of events
       #this is a list of the number of events of each type
       nevents = np.random.poisson(lam = mevs)
       for i, thisn in enumerate(nevents): 
          if thisn > 0:
             S, P1, P2, P3, B = events.get(i)(int(thisn), S, P1, P2, P3, B, t, params)
    return S, P1, P2, P3, B, t+dt

#one realisation of the model to get an incubation period
def onerun(thisN, params):
    
    # thisN = initial number of spores to be phagocytosed and that will lead to rupture
    S, P1, P2, P3, B = thisN, 0, 0, 0, 0
    t = 0.0
    Bmax = int(np.log(10**-5)/np.log(params['mu']/params['lamb']) + 1)
    popsum = S + P1 + P2 + P3 + B
    
    while B<Bmax and popsum>0:
       S, P1, P2, P3, B, t = onestep(t, S, P1, P2, P3, B, params)
       popsum = int(S) + int(P1) + int(P2) + int(P3) + int(B)
    
    if popsum == 0: #if infection goes extinct, incubation time is infinite
        return float('inf')
    
    else: #if number of bacteria reaches threshold Bmax, run deterministic part of the model
        
        initial_cond = [S, P1, P2, P3, B] #initial conditions for deterministic part of the model
        
        days = 0
        det_times = [0]
        #initiate symtoms_index to be equal to len(det_times) so the while loop starts
        symptoms_index = 1
        while symptoms_index == len(det_times):
            #start with giving it 5 days
            days+=5
            #times up to some number of days, with 12 time points every hour
            det_times = np.linspace(0, 24*days, 12*24*days+1)
            #get the output of the determninistic model
            B = get_ODE_solution(initial_cond, det_times, params)
            threshold_B = 10**10
            symptoms_index = np.searchsorted(B, threshold_B)
            #if symptoms_index == len(det_times) then the while loop will do it again with a bigger number of days
        #time between reaching Bmax and symtoms onset
        TNS = det_times[symptoms_index]
        infection_time = t + TNS
    return infection_time
    
#multiple realisations of the model to get a sample of finite incubation periods
def simulation(params, Nrealisations, inhaled_dose, sverdlovsk_dose = True):
    
    #sigma = lamb-mu
    #p = mu/(mu+lamb)
    params['lamb']  = params['sigma']*(1-params['p'])/(1-2*params['p'])
    params['mu'] = params['sigma']*params['p']/(1-2*params['p'])
    
    #phi = phi_hat*q is chosen to give an exponential dose-response curve with an ID50 of 8600
    phi = params['r']*(1+params['p']/(params['R']*(1-2*params['p'])))
    
    infection_times = []
    for _ in range(Nrealisations):
        if sverdlovsk_dose:
            #dose is sampled from the dose distribution estimated
            #for the Sverdlovsk outbreak
            inhaled_dose = sample_dose()
        inf_time = float('inf')
        while inf_time == float('inf'):
            #inhaled_dose will automatically be truncated to an integer for the binomial sampling
            deposited_dose = np.random.binomial(inhaled_dose, phi)
            inf_time = onerun(deposited_dose, params)
        infection_times.append(inf_time)
    model_incubation_times = np.array(infection_times)/24
    counts_m, bins_m = np.histogram(model_incubation_times, bins = np.linspace(0, 41, 42))

    return model_incubation_times, np.cumsum(counts_m)/len(model_incubation_times)