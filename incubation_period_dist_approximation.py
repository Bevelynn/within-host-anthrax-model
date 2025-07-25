import numpy as np
import matplotlib.pyplot as plt
import math
import human_model as model

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'axes.labelsize': 'x-large'})
plt.rcParams.update({'axes.titlesize': 'large'})
plt.rcParams.update({'xtick.labelsize': 'large'})
plt.rcParams.update({'ytick.labelsize': 'large'})

params = {'r': 8.06*10**-5, 'R': 1.6, 'rho': 0.004, 'delta': 0.107, 'sigma': 0.217, 'p': 0.01, 'N': 10**10}

#returns an array of the prob of symptom onset before time t, for all t in times (in hours),
#given dose D
def F(D, times, params, conditional = False):
    
    #extracellular bacterial replication rate
    lamb = params['sigma']*(1-params['p'])/(1-2*params['p'])
    #extracellular bacterial death rate
    mu = params['sigma']*params['p']/(1-2*params['p'])

    #N is the threshold of number of bacteria assumed to correspond to symtoms onset
    T = np.log(params['N'])/(lamb-mu)
    #print('T =', T/24)
    
    F1 = []
    for t in times:
        #prob of symptom onset before time t, starting with one spore,
        #given infection establishment
        if t <= T:
            F1.append(0)
        else:
            F1.append(Fr(t-T, params['rho'], params['delta']))
    if conditional:
        return (1-np.exp(-D*params['r']*np.array(F1)))/(1-np.exp(-D*params['r']))
    return 1-np.exp(-D*params['r']*np.array(F1))

#prob that rupture event occurs before time t, starting with one spore,
#given infection establishment
def Fr(t, rho, delta):
    return 1 - np.exp(-rho*t)*(delta/(delta-rho))**3 - sum([1/math.factorial(k)*np.exp(-delta*t)*(delta*t)**k*(1-(delta/(delta-rho))**(3-k)) for k in range(3)])


infection_doses = [1, 10, 100, 1000, 10**4, 10**5]
Nrealisations = 400
max_time = 60 # max time (in days) to plot the cdf for
# plot the probability every 4 hours up to the specified max incubation time
plotting_times = np.linspace(0, max_time, max_time*6+1) # in days
cum_probs = []
approx_cdfs = []

for dose in infection_doses:
    
    print('running realisations for dose', dose)
    #obtaining the cumulative fractions of incubation times from simulations
    inf_times = model.simulation(params, Nrealisations, dose, sverdlovsk_dose = False)[0]
    counts_d, bins_d = np.histogram(inf_times, bins = plotting_times)
    probs = np.cumsum(counts_d)/Nrealisations
    cum_probs.append(probs)
    
    #obtaining the approximate cdf of the incubation time
    cdf = F(dose, plotting_times[1:]*24, params, conditional = True)
    approx_cdfs.append(cdf)

dose_titles = ['1', '10', '100', '$10^3$', '$10^4$', '$10^5$']
plt.subplots(2, 3, figsize = (20, 10))
plt.subplots_adjust(left = None, bottom = None, right = None, top = 1.24, wspace = 0.2, hspace = 0.2)
for i, dose in enumerate(infection_doses):
    plt.subplot(2, 3, i+1)
    plt.plot(plotting_times[1:], cum_probs[i], linewidth = 5, label='Stochastic simulations')
    plt.plot(plotting_times[1:], approx_cdfs[i], linewidth = 5, ls = '--', label = 'Approximation')
    plt.title('Dose = %s' %dose_titles[i])
    
    if i in [3, 4, 5]:
        plt.xlabel('Incubation period (days)')
    else:
        plt.xticks([])
        
    if i in [0, 3]:
        plt.ylabel('Cumulative probability')
    else:
        plt.yticks([])
        
    if i == 0:
        plt.legend(ncols = 2, loc = (0, 1.1))
        
plt.savefig('Fig13.png', bbox_inches = 'tight')



