import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas
import animal_model as model

import warnings
# Suppress FutureWarning messages
warnings.simplefilter(action ='ignore', category = FutureWarning)

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'axes.labelsize': 'x-large'})
plt.rcParams.update({'axes.titlesize': 'large'})
plt.rcParams.update({'xtick.labelsize': 'large'})
plt.rcParams.update({'ytick.labelsize': 'large'})

########################################################################
'''##### load and plot posterior distributions for each species #####'''
########################################################################

sorted_params_rabbit = np.loadtxt('sorted_posterior_rabbits.csv', delimiter = ',')
#the order of the parameters here is:
    #q, delta, lamb_LN, mu_LN, M, mLN, lambC (fixed), mC
sorted_params_gp = np.loadtxt('sorted_posterior_guinea_pigs.csv', delimiter = ',')
#the order of the parameters here is:
    #q, delta (fixed), lamb_LN (fixed), mu_LN, M, mLN (fixed), lambC, mC, beta, muT
    
rabbit_sample_size = sorted_params_rabbit.shape[0]
print('rabbit sample size =', rabbit_sample_size)
gp_sample_size = sorted_params_gp.shape[0]
print('guinea-pig sample size =', gp_sample_size)

numb = 10**6
priorq = np.linspace(-5, 0, numb)[:-1]
priordelta = np.linspace(-4, 0.5, numb)[:-1]
priorlamb_LN = np.linspace(-1, 0.5, numb)[:-1]
priormu_LN = np.linspace(-4, 0.5, numb)[:-1]
priorM = np.linspace(1, 6, numb)[:-1]
priormLN = np.linspace(-4, 0.5, numb)[:-1]
priorlambC = np.linspace(-1, 0.5, numb)[:-1]
priormC = np.linspace(-1, 1, numb)[:-1]
priorbeta = np.linspace(-7,0,numb)[:-1]
priormuT = np.linspace(-3,1,numb)[:-1]

labels = ['$\log_{10}q$', '$\log_{10}\delta  \ (h^{-1})$', '$\log_{10}\lambda_{LN} \ (h^{-1})$', '$\log_{10}\mu_{LN} \ (h^{-1})$', '$\log_{10}M$ (CFU)', '$\log_{10}m_{LN} \ (h^{-1})$', '$\log_{10}\lambda_{C} \ (h^{-1})$', '$\log_{10}m_C \ (h^{-1})$', r'$\log_{10}\beta \ ng (CFU \cdot h)^{-1}$', '$\log_{10}\mu_T \ (h^{-1})$']
priors = [priorq, priordelta, priorlamb_LN, priormu_LN, priorM, priormLN, priorlambC, priormC, priorbeta, priormuT]

y, x = 2, 5

plt.subplots(y, x, figsize = (24,6))
plt.subplots_adjust(left = None, bottom = None, right = None, top = 1.24, wspace = 0.7, hspace = 0.6)

for param_index in [8,9]:
    
    plt.subplot(y, x, param_index+1)
    sns.kdeplot(priors[param_index], fill = True, color = 'gray', alpha = 0.6, label = 'prior')
    sns.kdeplot(sorted_params_gp[:,param_index], fill = True, color = 'orange', alpha = 0.6, label = 'guinea pig posterior')
    plt.grid()
    plt.xlabel(labels[param_index])
    
for param_index in range(8):
    
    plt.subplot(y, x, param_index+1)
    sns.kdeplot(priors[param_index], fill = True, color = 'gray', alpha = 0.6, label = 'prior')
    
    if param_index == 6:
        plt.plot([np.log10(0.17), np.log10(0.17)], [0, 1.25], color = 'blue')
    else:
        sns.kdeplot(sorted_params_rabbit[:,param_index], fill = True, color = 'blue', alpha = 0.6, label = 'rabbit posterior')
    
    if param_index == 1:
        plt.plot([np.log10(0.02), np.log10(0.02)], [0, 9], color = 'orange')
    elif param_index == 2:
        plt.plot([np.log10(1.49), np.log10(1.49)], [0, 5], color = 'orange')
    elif param_index == 5:
        plt.plot([np.log10(0.24), np.log10(0.24)], [0, 1], color = 'orange')
    else:
        sns.kdeplot(sorted_params_gp[:,param_index], fill = True, color = 'orange', alpha = 0.6, label = 'guinea pig posterior')
    
    plt.grid()
    plt.xlabel(labels[param_index])
    if param_index == 0:
        plt.legend(ncols = 3, loc = (0, 1.1), fontsize = 24)
        
plt.savefig('Fig4.png', bbox_inches = 'tight')


##################################
'''##### load rabbit data #####'''
##################################

rabbit_data_TBLN = pandas.read_excel('rabbit_data.xlsx', sheet_name = 0, header = 0)
TBLN_plotting_times = rabbit_data_TBLN['times'].to_numpy()
TBLN_plotting_CFU = rabbit_data_TBLN['CFU observed'].to_numpy()
TBLN_plotting_times_greater = rabbit_data_TBLN['times lower bound'].to_numpy()
TBLN_plotting_CFU_greater = rabbit_data_TBLN['CFU lower bound'].to_numpy()

rabbit_data_blood = pandas.read_excel('rabbit_data.xlsx', sheet_name = 1, header = 0)
blood_plotting_times = rabbit_data_blood['times'].to_numpy()
blood_plotting_CFU = rabbit_data_blood['CFU observed'].to_numpy()
blood_plotting_times_greater = rabbit_data_blood['times lower bound'].to_numpy()
blood_plotting_CFU_greater = rabbit_data_blood['CFU lower bound'].to_numpy()

rabbit_dose_response_data = pandas.read_excel('rabbit_data.xlsx', sheet_name = 2, header = 0)
rabbit_data_doses = rabbit_dose_response_data['dose'].to_numpy()
rabbit_data_probs = rabbit_dose_response_data['died'].to_numpy()/rabbit_dose_response_data['total animals'].to_numpy()


####################################################################################
'''##### get best parameter set for rabbits and plot stochastic simulations #####'''
####################################################################################

params_rabbit = {'phi_hat': 0.092, 'rho': 0.0735, 'R': 1.6, 'KLN': 10**9, 'KC':10**11, 'beta':0, 'muT':0}

best_params_rabbit = params_rabbit.copy()
best_params_rabbit['q'], best_params_rabbit['delta'], best_params_rabbit['lambLN'], best_params_rabbit['muLN'], best_params_rabbit['M'], best_params_rabbit['mLN'], best_params_rabbit['lambC'], best_params_rabbit['mC'] = pow(10, sorted_params_rabbit[-1,:])
best_params_rabbit['M'] = int(best_params_rabbit['M'])

print('rabbit best fit parameters', best_params_rabbit)

#granularity of time points
gran = 2 #one point every half hour
tmax = 40
sim_times = np.linspace(0, tmax, gran*tmax+1)
time_points = len(sim_times)

Nrealisations = 10  # the number of realisations to do
Sims_Bln = np.zeros((Nrealisations, time_points))
Sims_Bc = np.zeros((Nrealisations, time_points))
for i in range(Nrealisations):

    #sample the initial number of CFU in the lungs that will be 
    #phagocytosed and lead to a rupture in the lymph nodes
    thisCFU = np.random.binomial(model.rabbit_data_dose, best_params_rabbit['phi_hat']*best_params_rabbit['q'])
    
    #do one realisation
    realisation = model.onerun(thisCFU, sim_times, best_params_rabbit)
     
    #store the timcourse for each variable in the corresponding array
    Sims_Bln[i,] = realisation[1]
    Sims_Bc[i,] = realisation[2]
    
    print('finished best realisation %s' %i)
    
    
plt.subplots(1, 2, figsize = (24, 6))
plt.subplots_adjust(left = None, bottom = None, right = None, top = 1.24, wspace = 0.2, hspace = 0.2)

plt.subplot(1, 2, 1)

for i in range(len(Sims_Bln)):
    plt.plot(sim_times, Sims_Bln[i, :])
plt.scatter(TBLN_plotting_times, TBLN_plotting_CFU, color = 'black', zorder = 500)
plt.scatter(TBLN_plotting_times_greater, TBLN_plotting_CFU_greater, color = 'black', marker = 's', zorder = 500)

plt.xlabel('Time (hours)')
plt.ylabel('Bacterial CFU')
plt.ylim(1,10**10)
plt.yscale('log')
plt.xticks(np.arange(0, 48, 12))
plt.title('CFU in lymph nodes')

plt.subplot(1, 2, 2)

for i in range(len(Sims_Bc)):
    plt.plot(sim_times, Sims_Bc[i, :])
plt.scatter(blood_plotting_times, blood_plotting_CFU, color = 'black', zorder = 500)
plt.scatter(blood_plotting_times_greater, blood_plotting_CFU_greater, color = 'black', marker = 's', zorder = 500)

plt.xlabel('Time (hours)')
plt.ylim(1,10**10)
plt.yscale('log')
plt.xticks(np.arange(0, 48, 12))
plt.title('CFU in blood')

plt.savefig('Fig2.png', bbox_inches = 'tight')


#############################################################
'''##### plot rabbit dose-response curve predictions #####'''
#############################################################

#doses to plot the predicted dose-response curve for
doses = np.logspace(2,7,50)
#empty array to add the predicted dose-response curves for all posterior parameter sets
dose_response_probs = np.empty((0,len(doses)))

for i in range(rabbit_sample_size):
    params_rabbit['q'], params_rabbit['delta'], params_rabbit['lambLN'], params_rabbit['muLN'], params_rabbit['M'], params_rabbit['mLN'], params_rabbit['lambC'], params_rabbit['mC'] = pow(10, sorted_params_rabbit[i, :])
    probs = model.get_dose_response(params_rabbit, doses)
    dose_response_probs = np.vstack((dose_response_probs, probs))

probs_best = dose_response_probs[-1,:]

#Calculate the pointwise median and credible intervals of the probabilities for each dose
medians = [np.median(dose_response_probs[:,i]) for i in range(len(doses))]
low_CIs = [np.percentile(dose_response_probs[:,i],2.5) for i in range(len(doses))]
high_CIs = [np.percentile(dose_response_probs[:,i],97.5) for i in range(len(doses))]

plt.figure(figsize = (8, 7))
plt.scatter(rabbit_data_doses, rabbit_data_probs)
plt.plot(doses, medians, label='pointwise median')
plt.fill_between(doses, high_CIs, low_CIs, alpha = 0.3, label = '95$\%$ CI', color = 'steelblue')
plt.xticks([10, 100, 1000, 10**4, 10**5, 10**6, 10**7], labels = ['$10$','$10^2$','$10^3$','$10^4$','$10^5$','$10^6$','$10^7$'])
plt.xscale('log')
plt.xlabel('Inhaled dose of spores')
plt.ylabel('Probability of infection')
plt.legend(fontsize = 18)
plt.savefig('Fig3.png', bbox_inches='tight')


######################################
'''##### load guinea pig data #####'''
######################################

gp_data_CFU = pandas.read_excel('guinea_pig_data.xlsx', sheet_name = 0, header = 0)
CFU_plotting_times = gp_data_CFU['times'].to_numpy()
CFU_plotting = gp_data_CFU['CFU observed'].to_numpy()
CFU_plotting_times_greater = gp_data_CFU['times lower bound'].to_numpy()
CFU_plotting_greater = gp_data_CFU['CFU lower bound'].to_numpy()
CFU_plotting_times_less = gp_data_CFU['times upper bound'].to_numpy()
CFU_plotting_less = gp_data_CFU['CFU upper bound'].to_numpy()

gp_data_PA = pandas.read_excel('guinea_pig_data.xlsx', sheet_name = 1, header = 0)
PA_plotting_times = gp_data_PA['times'].to_numpy()
PA_plotting = gp_data_PA['PA observed'].to_numpy()
PA_plotting_times_greater = gp_data_PA['times lower bound'].to_numpy()
PA_plotting_greater = gp_data_PA['PA lower bound'].to_numpy()
PA_plotting_times_less = gp_data_PA['times upper bound'].to_numpy()
PA_plotting_less = gp_data_PA['PA upper bound'].to_numpy()

gp_dose_response_data = pandas.read_excel('guinea_pig_data.xlsx', sheet_name = 2, header = 0)
gp_data_doses = gp_dose_response_data['dose'].to_numpy()
gp_data_probs = gp_dose_response_data['died'].to_numpy()/gp_dose_response_data['total animals'].to_numpy()


########################################################################################
'''##### get best parameter set for guinea pigs and plot stochastic simulations #####'''
########################################################################################

params_gp = {'phi_hat': 0.3, 'rho': 0.0735, 'R': 1.6, 'KLN': 10**9, 'KC':10**11}

best_params_gp = params_gp.copy()
best_params_gp['q'], best_params_gp['delta'], best_params_gp['lambLN'], best_params_gp['muLN'], best_params_gp['M'], best_params_gp['mLN'], best_params_gp['lambC'], best_params_gp['mC'], best_params_gp['beta'], best_params_gp['muT'] = pow(10, sorted_params_gp[-1,:])
best_params_gp['M'] = int(best_params_gp['M'])

print('guinea pig best fit parameters', best_params_gp)

#granularity of time points
gran = 2 #one point every half hour
tmax = 75
sim_times = np.linspace(0, tmax, gran*tmax+1)
time_points = len(sim_times)

Nrealisations = 10 # the number of realisations to do
Sims_Bln = np.zeros((Nrealisations, time_points))
Sims_Bc = np.zeros((Nrealisations, time_points))
Sims_PAln = np.zeros((Nrealisations, time_points))
Sims_PAc = np.zeros((Nrealisations, time_points))
for i in range(Nrealisations):

    #sample the initial number of CFU in the lungs that will be 
    #phagocytosed and lead to a rupture in the lymph nodes
    thisCFU = np.random.binomial(model.gp_data_dose, best_params_gp['phi_hat']*best_params_gp['q'])
    
    #do one realisation
    realisation = model.onerun(thisCFU, sim_times, best_params_gp)
     
    #store the timcourse for each variable in the corresponding array
    Sims_Bln[i,] = realisation[1]
    Sims_Bc[i,] = realisation[2]
    Sims_PAln[i,] = realisation[3]
    Sims_PAc[i,] = realisation[4]
    
    print('finished best realisation %s' %i)
    
    
plt.subplots(2, 2, figsize = (24, 12))
plt.subplots_adjust(left = None, bottom = None, right = None, top = 1.24, wspace = 0.2, hspace = 0.2)

plt.subplot(2, 2, 1)

for i in range(len(Sims_Bln)):
    plt.plot(sim_times, Sims_Bln[i, :])

plt.ylabel('Bacterial CFU')
plt.ylim(1, 10**12)
plt.yscale('log')
plt.xticks(np.arange(0, 75, 12))
plt.title('Lymph nodes')

plt.subplot(2, 2, 2)

for i in range(len(Sims_Bc)):
    plt.plot(sim_times, Sims_Bc[i, :])
plt.scatter(CFU_plotting_times, CFU_plotting, color = 'black', zorder = 500)
plt.scatter(CFU_plotting_times_greater, CFU_plotting_greater, color = 'black', marker = 's', zorder = 500)
plt.scatter(CFU_plotting_times_less, CFU_plotting_less, color = 'black', marker = 'v', zorder = 500)
plt.plot(sim_times,[CFU_plotting_less]*len(sim_times), '--', color = 'black', alpha = 0.5, label = 'LLOD')

plt.ylim(1, 10**12)
plt.yscale('log')
plt.xticks(np.arange(0, 75, 12))
plt.title('Blood')

plt.subplot(2, 2, 3)

for i in range(len(Sims_PAln)):
    plt.plot(sim_times, Sims_PAln[i, :])

plt.xlabel('Time (hours)')
plt.ylabel('Amount of PA (ng)')
plt.ylim(1, 10**12)
plt.yscale('log')
plt.xticks(np.arange(0, 75, 12))


plt.subplot(2, 2, 4)

for i in range(len(Sims_PAc)):
    plt.plot(sim_times, Sims_PAc[i, :])
plt.scatter(PA_plotting_times, PA_plotting, color = 'black', zorder = 500)
plt.scatter(PA_plotting_times_greater, PA_plotting_greater, color = 'black', marker = 's', zorder = 500)
plt.scatter(PA_plotting_times_less, PA_plotting_less, color = 'black', marker = 'v', zorder = 500)
plt.plot(sim_times,[PA_plotting_less[0]]*len(sim_times), '--', color = 'black', alpha = 0.5, label = 'LLOD')

plt.xlabel('Time (hours)')
plt.ylim(1, 10**12)
plt.yscale('log')
plt.xticks(np.arange(0, 75, 12))

plt.savefig('Fig5.png', bbox_inches = 'tight')


#################################################################
'''##### plot guinea pig dose-response curve predictions #####'''
#################################################################

#doses to plot the predicted dose-response curve for
doses = np.logspace(2, 7, 50)
#empty array to add the predicted dose-response curves for all posterior parameter sets
dose_response_probs = np.empty((0, len(doses)))

for i in range(gp_sample_size):
    params_gp['q'], params_gp['delta'], params_gp['lambLN'], params_gp['muLN'], params_gp['M'], params_gp['mLN'], params_gp['lambC'], params_gp['mC'], params_gp['beta'], params_gp['muT'] = pow(10, sorted_params_gp[i, :])
    probs = model.get_dose_response(params_gp, doses)
    dose_response_probs = np.vstack((dose_response_probs, probs))

probs_best = dose_response_probs[-1,:]

#Calculate the pointwise median and credible intervals of the probabilities for each dose
medians = [np.median(dose_response_probs[:,i]) for i in range(len(doses))]
low_CIs = [np.percentile(dose_response_probs[:,i],2.5) for i in range(len(doses))]
high_CIs = [np.percentile(dose_response_probs[:,i],97.5) for i in range(len(doses))]


plt.figure(figsize = (8, 7))
plt.scatter(gp_data_doses, gp_data_probs)
plt.plot(doses, medians, label = 'pointwise median')
plt.fill_between(doses, high_CIs, low_CIs, alpha=0.3, label = '95$\%$ CI', color = 'steelblue')
plt.xticks([10, 100, 1000, 10**4, 10**5, 10**6, 10**7], labels = ['$10$','$10^2$','$10^3$','$10^4$','$10^5$','$10^6$','$10^7$'])
plt.xscale('log')
plt.xlabel('Inhaled dose of spores')
plt.ylabel('Probability of infection')
plt.legend(fontsize = 18)
plt.savefig('Fig6.png', bbox_inches = 'tight')
