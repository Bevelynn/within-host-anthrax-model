import numpy as np
from scipy.integrate import solve_ivp

#timestep for tau leaping
tau_dt = 0.01

#mean inhaled dose of spores to which animals were exposed in the experiments
rabbit_data_dose = 4.428*10**7
gp_data_dose = 2*10**7

#################################################
'''functions to simulate the model'''
#################################################


#differential equation for PA levels, given a fixed number of bacteria in the lymph nodes and circulation
def diff_eqns_PA(t, x, BLN, BC, params):
    
    PALN, PAC = x
    
    dPALN = params['beta']*BLN - (params['muT'] + params['mLN'])*PALN
    dPAC = params['mLN']*PALN + params['beta']*BC - (params['muT'] + params['mC'])*PAC
    
    return np.hstack((dPALN, dPAC))

#solve the PA eqns to get the amount of PA at the next timestep
def solution_PA(BLNold, BCold, PALNold, PACold, timestep, params):
    
    sol = solve_ivp(fun = lambda t, y: diff_eqns_PA(t, y, BLNold, BCold, params), t_span=[0,timestep], y0=[PALNold, PACold], method='BDF', t_eval=[timestep])
    
    PA_LN = sol.y[0,:]
    PA_C = sol.y[1,:]

    return PA_LN[0], PA_C[0]

#differential equations for each variable of the model
def derivative(t, x, params):
    #variables are:
    #spores in airways that will lead to infected cells that rupture
    #infected cells in stage 1-3
    #CFU in lymph nodes
    #CFU in blood
    #PA in lymph nodes
    #PA in blood
    
    S, P1, P2, P3, BLN, BC, PALN, PAC = x
    
    dS = -params['rho']*S
    dP1 = params['rho']*S - params['delta']*P1
    dP2 = params['delta']*(P1-P2)
    dP3 = params['delta']*(P2-P3)
    dBLN = params['delta']*params['R']*P3 + params['lambLN']*BLN*(1-BLN/params['KLN']) - (params['muLN'] + params['mLN'])*BLN
    dBC = params['mLN']*BLN - params['mC']*BC + params['lambC']*BC*(1-BC/params['KC'])
    dPALN = params['beta']*BLN - (params['muT'] + params['mLN'])*PALN
    dPAC = params['mLN']*PALN + params['beta']*BC - (params['muT'] + params['mC'])*PAC
    
    return np.hstack((dS, dP1, dP2, dP3, dBLN, dBC, dPALN, dPAC))

#solve the model eqns
def get_ODE_solution(initial_cond, det_times, params):

    sol = solve_ivp(fun = lambda t, x: derivative(t, x, params), t_span=[det_times[0], det_times[-1]], y0=initial_cond, method='BDF', t_eval=det_times)
    
    B_LN = sol.y[4,:]
    B_C = sol.y[5,:]
    PA_LN = sol.y[6,:]
    PA_C = sol.y[7,:]

    return B_LN, B_C, PA_LN, PA_C
   
#functions for the different events that can happen in the stochastic process
def M_eat_S_n(n, S, P1, P2, P3, BLN, BC, params):
    """ macrophages phagocytose n spores """
    n = min(n,S)
    return S-n, P1+n, P2, P3, BLN, BC

def M_1_2_n(n, S, P1, P2, P3, BLN, BC, params):
    """ n macrophages go from stage 1 to stage 2 """
    n = min(n,P1)
    return S, P1-n, P2+n, P3, BLN, BC

def M_2_3_n(n, S, P1, P2, P3, BLN, BC, params):
    """ n macrophages go from stage 2 to stage 3 """
    n = min(n,P2)
    return S, P1, P2-n, P3+n, BLN, BC

def M_rupt_n(n, S, P1, P2, P3, BLN, BC, params):
   """ n host macrophages rupture and release bacteria into the lymph nodes"""
   for _ in range(min(n,P3)):
       ruptsize = np.random.geometric(1/params['R'])
       P3-=1
       BLN+=ruptsize
   return S, P1, P2, P3, BLN, BC

def B_birth_LN_n(n, S, P1, P2, P3, BLN, BC, params):
    """ n extracellular bacteria in the lymph nodes replicate """
    return S, P1, P2, P3, BLN+n, BC

def B_death_LN_n(n, S, P1, P2, P3, BLN, BC, params):
    """ n extracellular bacteria in the lymph nodes die """
    return S, P1, P2, P3, max(0,BLN-n), BC

def B_migrate_n(n, Se, P1, P2, P3, BLN, BC, params):
    """ n extracellular bacteria in the lymph nodes migrate to the blood """
    n = min(n, BLN)
    return Se, P1, P2, P3, BLN-n, BC+n

events={0:M_eat_S_n, 1:M_1_2_n, 2:M_2_3_n, 3:M_rupt_n, 4:B_birth_LN_n, 5:B_death_LN_n, 6:B_migrate_n}

#rate of each event
def makerates(S, P1, P2, P3, BLN, params):
   
   myrates = []
   
   myrates.append(params['rho']*S)                                              # phagocytosis of spores
   myrates.append(params['delta']*P1)                                           # next stage for infected cells
   myrates.append(params['delta']*P2)                                           # next stage for infected cells
   myrates.append(params['delta']*P3)                                           # next stage for infected cells
   myrates.append(params['lambLN']*BLN)                                         # replication of extrac bacteria in LN
   myrates.append(params['lambLN']*BLN**2/params['KLN'] + params['muLN']*BLN)   # death of extrac bacteria in LN
   if BLN > params['M']:
       myrates.append(params['mLN']*BLN)                                        # migration of extrac bacteria into C
   else:
       myrates.append(0)
   
   return myrates

#store the values of each variable in lists
def writestuff(t_index, BLN, BC, PALN, PAC, timecourses):
    
    for idx, val in enumerate([t_index, BLN, BC, PALN, PAC]):
        timecourses[idx].append(val)
    
    return timecourses
   
#one step of the Gillespie or tau-leaping algorithm
def onestep(t, S, P1, P2, P3, BLN, BC, PALN, PAC, timecourses, sim_times, params):
    
    rates = makerates(S, P1, P2, P3, BLN, params)
    BLNold, BCold, PALNold, PACold = BLN, BC, PALN, PAC
    totalrate = sum(rates)
    #if the average time to do one event is very small (e.g less than 0.001) then do tau leaping
    if 1/totalrate<0.001:
        useGillespie = False
    else:
        useGillespie = True
    if useGillespie:
       urv = np.random.uniform()
       i = np.searchsorted(np.cumsum(rates),urv*totalrate)
       S, P1, P2, P3, BLN, BC = events.get(i)(1, S, P1, P2, P3, BLN, BC, params)
       dt = -np.log(np.random.uniform())/totalrate
    else:
       dt = tau_dt
       mevs = [dt*rate for rate in rates]  # mean numbers of events
       #choose number of events from poisson dist with lambda=mean number of events
       #this is a list of the number of events of each reaction type
       nevents = np.random.poisson(lam = mevs)
       for i,thisn in enumerate(nevents): 
          if thisn > 0:
             S, P1, P2, P3, BLN, BC = events.get(i)(int(thisn), S, P1, P2, P3, BLN, BC, params)
    
    PALN, PAC = solution_PA(BLNold, BCold, PALNold, PACold, dt, params)
    
    #only save the model variables at specific times specified by sim_times
    #so that the time series of each run is the same length and we can match
    #them up easily when plotting
    i = 1
    prev_i = timecourses[0][-1]
    #if the next time in sim_times is before the time of the event that has
    #just happened, the process will still be at the previous state at
    #that time, so we need to write the previous state to the time series lists
    while (prev_i+i)<len(sim_times) and t+dt>sim_times[prev_i+i]:
       writestuff(prev_i+i, BLNold, BCold, PALNold, PACold, timecourses)
       i+=1
       
    return S, P1, P2, P3, BLN, BC, PALN, PAC, t+dt, timecourses

#one full realisation of the model
def onerun(thisN, sim_times, params):
    
    #thisN is the initial number of spores to be phagocytosed and that will lead to rupture
    S, P1, P2, P3, BLN, BC, PALN, PAC = thisN, 0, 0, 0, 0, 0, 0, 0
    t = 0.0
    timecourses = [[0], [0], [0], [0], [0]]
    popsum = S + P1 + P2 + P3 + BLN
    
    #run stochastic simulations until either
    #bacteria enters circulation or all spore/bacteria populations are extinct
    while BC < 1 and popsum > 0:
       S, P1, P2, P3, BLN, BC, PALN, PAC, t, timecourses = onestep(t, S, P1, P2, P3, BLN, BC, PALN, PAC, timecourses, sim_times, params)
       popsum = int(S) + int(P1) + int(P2) + int(P3) + int(BLN)
    
    #prev_i is the index of the last time point that was written to the timecourses
    prev_i = timecourses[0][-1]
    
    if popsum == 0: #if extinction occurs, we keep increasing the index until we have filled the timecourse for all times
        i=1
        while (prev_i+i)<len(sim_times):
           #if the time series has not been filled yet, fill in the current
           #model state (i.e. extinction) at the next time in sim_times
           #note this ignores natural PA degradation after all other populations are zero
           writestuff(prev_i+i, BLN, BC, PALN, PAC, timecourses)
           i+=1
    
    else: #if BC>=1 (bacteria has entered circulation), we simulate the deterministic model for all later times
        if (prev_i+1)<len(sim_times):
            #if the time series has not been filled yet, simulate the deterministic 
            #model starting from time t and evaluated at all later times in sim_times
            det_times = np.hstack((np.array(t),sim_times[prev_i+1:]))
            initial_cond = [S, P1, P2, P3, BLN, BC, PALN, PAC]
            det_sim = get_ODE_solution(initial_cond, det_times, params)
            bacLN, bacC, paLN, paC = det_sim
            #write the values of BLN, BC, PALN, and PAC to the timecourses for all times in det_times that are also in sim_times
            #i.e. every time apart from the first one of det_times
            for j in range(1, len(det_times)):
                writestuff(prev_i+j, bacLN[j], bacC[j], paLN[j], paC[j], timecourses)
    
    return timecourses


################################################
'''funtion to compare with dose-response data'''
#################################################

#probability of infection for a list of inhaled doses
def get_dose_response(params, doses):

    p = params['muLN']/(params['lambLN'] + params['muLN'])
    r = params['phi_hat']*params['q']*params['R']*(1-2*p)/(params['R']*(1-2*p)+p)
    prob_infs = 1-np.exp(-r*doses)
    
    return prob_infs