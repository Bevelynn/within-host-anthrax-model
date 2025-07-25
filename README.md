## Description of .xlsx and .csv files: 
 
### "rabbit_data.xlsx":
Required for the file "animal_calibration_visualiser.py".  
Contains experimental data used to calibrate the model for rabbits.  
Sheet 1: Data from Table 4 of CFU in TBLN.  
Sheet 2: Data from Table 4 of CFU in blood.  
Sheet 3: Rabbit dose-response data from Ref. \[22\].

### "guinea_pig_data.xlsx":
Required for the file "animal_calibration_visualiser.py".  
Contains experimental data used to calibrate the model for guinea pigs.  
Sheet 1: Data from Table 5 of CFU in blood.  
Sheet 2: Data from Table 5 of PA in blood.  
Sheet 3: Guinea-pig dose-response data from Ref. \[22\].

### "sorted_posterior_rabbits.csv":
Required for the file "animal_calibration_visualiser.py".  
Contains lists of parameter values from the posterior distribution for rabbits, ordered from largest to smallest distance.  
Column 1: Values of $\log_{10}q$.  
Column 2: Values of $\log_{10}\delta$.    
Column 3: Values of $\log_{10}\lambda_{LN}$.    
Column 4: Values of $\log_{10}\mu_{LN}$.    
Column 5: Values of $\log_{10}M$.    
Column 6: Values of $\log_{10}m_{LN}$.  
Column 7: Values of $\log_{10}\lambda_C$ (fixed for rabbits).  
Column 8: Values of $\log_{10}m_C$.  

### "sorted_posterior_guinea_pigs.csv":
Required for the file "animal_calibration_visualiser.py".  
Contains lists of parameter values from the posterior distribution for guinea pigs, ordered from largest to smallest distance.    
Column 1: Values of $\log_{10}q$.  
Column 2: Values of $\log_{10}\delta$ (fixed for guinea pigs).  
Column 3: Values of $\log_{10}\lambda_{LN}$ (fixed for guinea pigs).  
Column 4: Values of $\log_{10}\mu_{LN}$.  
Column 5: Values of $\log_{10}M$.    
Column 6: Values of $\log_{10}m_{LN}$ (fixed for guinea pigs).  
Column 7: Values of $\log_{10}\lambda_C$.  
Column 8: Values of $\log_{10}m_C$.  
Column 9: Values of $\log_{10}\beta$.  
Column 10: Values of $\log_{10}\mu_T$.  

### "Sverdlovsk_incubation_periods.csv":
Required for the files "Sverdlovsk_data_plots.py" and "human_calibration_visualiser.py".  
Contains a list of the number of days to symptom onset for 30 individuals from the Sverdlovsk outbreak.

### "Sverdlovsk_dose_distribution.csv":
Required for the file "human_model.py".  
Distribution of doses of spores to which individuals are assumed to have been exposed in the Sverdlovsk outbreak, estimated by Wilkening for their 'Model D' in Ref. \[31\].

### "sorted_posterior_human.csv":  
Required for the file "human_calibration_visualiser.py".  
Contains lists of parameter values from the posterior distribution for humans, ordered from largest to smallest distance.  
Column 1: Values of $\log_{10}\rho$.    
Column 2: Values of $\log_{10}\delta$.    
Column 3: Values of $\log_{10}\sigma$.    
Column 4: Values of $\log_{10}p$.    

          
## Description of .py files

### "animal_model.py":
Contains functions to simulate the within-host model for animals and a function to obtain the dose-response curve.

### "animal_calibration_visualiser.py":
Plots the posterior distributions in "sorted_posterior_rabbits.csv" and "sorted_posterior_guinea_pigs.csv", creating Figure 4.  
Uses "animal_model.py" to simulate stochastic realisations of the within-host model, using the parameter set with the smallest distance from the data for each species. Ten stochastic realisations are plotted for each species to compare with the data used, creating Figures 2 and 5.  
Obtains the pointwise median and 95% credible intervals of the dose-response curve predictions using all posterior parameter sets for each species to compare with the dose-response data used. These are used to produce Figures 3 and 6.

### "Sverdlovsk_data_plots.py":  
Plots the Sverdlovsk incubation period data set to produce Figure 8.

### "human_model.py":
Contains functions to simulate the within-host model for humans to obtain samples of incubation periods.

### "human_calibration_visualiser.py":
Uses "human_model.py" to simulate incubation period distributions for each parameter set in the posterior sample, "sorted_posterior_human.csv", and obtains the pointwise median and 95% credible intervals of these distributions. These are plotted with the Sverdlovsk incubation period data set to produce Figure 9.  
Plots the marginal posterior distributions from "sorted_posterior_human.csv", creating Figure 10.  
Transforms the parameter values in "sorted_posterior_human.csv" to obtain the posterior distributions of $\lambda$, $\mu$, and $\phi$, creating Figure 11.  
Produces scatter plots of pairs of parameter values from "sorted_posterior_human.csv" to indicate pairwise correlations, creating Figure 12.

### "incubation_period_dist_approximation.py":
Uses "human_model.py" to simulate incubation period distributions for a range of inhaled doses and compares these with the calculation of the approximate cumulative distribution function of the incubation period, creating Figure 13.
