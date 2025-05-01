import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas
import human_model as model

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'axes.labelsize': 'x-large'})
plt.rcParams.update({'axes.titlesize': 'large'})
plt.rcParams.update({'xtick.labelsize': 'large'})
plt.rcParams.update({'ytick.labelsize': 'large'})

##### Load incubation times data #####

incubation_times = pandas.read_csv('Sverdlovsk_incubation_periods_data.csv', header=0)['incubation time (days)'].to_numpy()
Nrealisations = len(incubation_times)

#obtaining the cumulative daily fractions of incubation times for days 0 to 40
counts_d, bins_d = np.histogram(incubation_times, bins = np.linspace(0,41,42))
data = np.cumsum(counts_d)/len(incubation_times)


##### Load posterior distribution #####

sorted_params = np.loadtxt('sorted_posterior_human.csv', delimiter = ',')
#the order of the parameters here is:
    #rho, delta, sigma, p
    
#number of parameter sets in the posterior sample
sample_size = sorted_params.shape[0]
#number of parameters estimated
nump = sorted_params.shape[1]

print('sample size =', sample_size)


##### Simulate an incubation period distribution for each accepted parameter set #####

params = {'r': 8.06*10**-5, 'R': 1.6}

time_points = 41
Sims = np.zeros((sample_size, time_points))
for i in range(sample_size):
   params['rho'], params['delta'], params['sigma'], params['p'] = pow(10, sorted_params[i,:])
   sim = model.simulation(params, Nrealisations, inhaled_dose = 'none')
   Sims[i,] = sim[1]
   print('finished realisations %s' %i)
   

##### Plot median and 95% credible interval of the simulation results vs. data used for the calibration #####

#calculating pointwise median and credible interval of the incubation period distributions
medians = [np.median(Sims[:,i]) for i in range(time_points)]
low_CIs = [np.percentile(Sims[:,i],2.5) for i in range(time_points)]
high_CIs = [np.percentile(Sims[:,i],97.5) for i in range(time_points)]

plt.figure(figsize = (10, 6))
plt.plot(np.linspace(0, 40, 41), medians, label = 'Pointwise median', color = 'blue')
plt.fill_between(np.linspace(0, 40, 41), high_CIs, low_CIs, alpha = 0.3, label = '95$\%$ CI', color = '#1f77b4')
plt.scatter(np.linspace(0, 40, 41), data, color = 'blue', label = 'Sverdlovsk data')
plt.legend(loc = 4, fontsize = 16)
plt.xlabel('Incubation period (days)')
plt.ylabel('Cumulative probability')
plt.savefig('Fig9.png', bbox_inches = 'tight')


##### Plot marginal prior and posterior densities #####

numb = 10**6

#Upper and lower bounds of the uniform priors for each of the parameters
#log10(rho), log10(delta), log10(sigma), log10(p)
lower_bounds = [-3, -3, -2, -4]
upper_bounds = [1, 1, np.log10(2), np.log10(0.49)]

priorrho = np.linspace(lower_bounds[0], upper_bounds[0], numb)[:-1]
priordelta = np.linspace(lower_bounds[1], upper_bounds[1], numb)[:-1]
priorsigma = np.linspace(lower_bounds[2], upper_bounds[2], numb)[:-1]
priorp = np.linspace(lower_bounds[3], upper_bounds[3], numb)[:-1]

labels = [r'$\log_{10}\rho \ (h^{-1})$','$\log_{10}\delta \ (h^{-1})$', '$\log_{10}\sigma \ (h^{-1})$', r'$\log_{10}p$']
priors = [priorrho, priordelta, priorsigma, priorp]

y, x = 1, 4

plt.subplots(y, x, figsize = (24,6))
plt.subplots_adjust(left = None, bottom = None, right = None, top = 1.24, wspace = 0.4, hspace = 0.4)

for param_index in range(nump):
    
    plt.subplot(y,x,param_index+1)
    sns.kdeplot(priors[param_index], fill = True, color = 'gray', alpha = 0.6)
    sns.kdeplot(sorted_params[:,param_index], fill = True, color = 'green', alpha = 0.6)
    plt.ylabel('Density', fontsize = 18)
    plt.grid()
    plt.xlabel(labels[param_index])

plt.tight_layout(w_pad = 1)
plt.savefig('Fig10.png', bbox_inches = 'tight')


###### Transforming to lambda, mu, and phi posteriors ######

sigma = pow(10, sorted_params[:,2])
p = pow(10, sorted_params[:,3])

mu = np.log10(sigma*p/(1-2*p))
lamb  = np.log10(sigma*(1-p)/(1-2*p))

phi = np.log10(params['r']*(1+p/(params['R']*(1-2*p))))

plt.subplots(1, 3, figsize = (18, 6))
plt.subplots_adjust(left = None, bottom = None, right = None, top = 1.24, wspace = 0.4, hspace = 0.4)

plt.subplot(131)
sns.kdeplot(lamb, fill = True, color = 'green', alpha = 0.6)
plt.ylabel('Density', fontsize = 18)
plt.grid()
plt.xlabel('$\log_{10}\lambda \ (h^{-1})$')

plt.subplot(132)
sns.kdeplot(mu, fill = True, color = 'green', alpha = 0.6)
plt.ylabel('Density', fontsize = 18)
plt.grid()
plt.xlabel('$\log_{10}\mu \ (h^{-1})$')

plt.subplot(133)
sns.kdeplot(phi, fill = True, color = 'green', alpha = 0.6)
plt.ylabel('Density', fontsize = 18)
plt.grid()
plt.xlabel('$\log_{10}\phi$')

plt.tight_layout(w_pad = 1)
plt.savefig('Fig11.png', bbox_inches = 'tight')


##### Scatterplots to indicate pairwise correlations #####

y, x = 2, 3

plt.subplots(y, x, figsize = (18, 12))
plt.subplots_adjust(left = None, bottom = None, right = None, top = 1.24, wspace = 0.6, hspace = 1)
taken = 0
for j in range(nump, 1, -1):
    for i in range(1, j):
        plt.subplot(y, x, i+taken)
        plt.scatter(sorted_params[:,nump-j], sorted_params[:,i + (nump-j)], color = 'green')
        plt.xlabel(labels[nump-j])
        plt.ylabel(labels[i + (nump-j)])
    taken+=j-1
plt.tight_layout(w_pad = 2)
plt.savefig('Fig12.png', bbox_inches = 'tight')
