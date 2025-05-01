import numpy as np
import matplotlib.pyplot as plt
import pandas

##############################
'''Plotting Sverdlovsk data'''
##############################

incubation_times = pandas.read_csv('Sverdlovsk_incubation_periods_data.csv', header = 0)['incubation time (days)'].to_numpy()

counts_d, bins_d = np.histogram(incubation_times, bins = np.linspace(0, 41, 42))


plt.subplots(1, 2, figsize = (15, 5))
plt.subplots_adjust(left = None, bottom = None, right = None, top = 1.24, wspace = 0.4, hspace = 0.4)

plt.subplot(121)
plt.hist(incubation_times, bins = np.linspace(0, 41, 42), align = 'left', density = bool, alpha = 0.5, edgecolor = 'blue', label = 'Sverdlovsk data')
plt.legend()
plt.xlabel('Incubation period (days)')
plt.ylabel('Density')

plt.subplot(122)
plt.bar(np.linspace(0, 41, 42)[:-1], np.cumsum(counts_d)/len(incubation_times), width = 1, alpha = 0.5, edgecolor = 'blue', label = 'Sverdlovsk data')
plt.xlabel('Incubation period (days)')
plt.ylabel('Cumulative probability')

plt.savefig('Fig8.png', bbox_inches = 'tight')