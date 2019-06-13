import numpy as np 
from src.da_methods import *
from src.pkg import *
from src.util import *
import matplotlib.pyplot as plt# Configure ES-MDA
nreaz = 100          # number of realizations
infl = 1.02          # inflation coefficient for state vector
niter = 4          # number of iterations
alpha = np.array([4,4,4,4])        # inflation coefficient for each iteration, SUM(1/alpha(i))=1 where i is the number of iteration; the coefficients could be the same or in decreasing order
da_para = ['flux', 'thermal conductivity','porosity']          # ['permeability','thermal conductivity','porosity','flux'],  'permeability' and 'flux' cannot be chosen together
da_timestep = 900 # unit: sec, data assimilation timestep

# Default mean, standard deviation, upper and lower boundaries for each variable
logperm_mean = -11          # unit: log(m2), mean of prior log(permeability)
logperm_sd = 2          # unit: log(m2), S.D. of prior log(permeability)
perm_low_bound = 1e-13          # unit: m2, lower bound of permeability
perm_up_bound = 1e-9          # unit: m2, upper bound of permeability

flux_mean = 0          # unit: m/day, mean flux
flux_sd = 0.5          # unit: m/day, standard deviation of flux
flux_low_bound = -5   # unit: m/day, lower bound of hydraulic gradient       
flux_up_bound = 5     # unit: m/day, upper bound of hydraulic gradient    

th_cond_mean = 1.5          # unit: W/(mK)-1, mean of prior thermal conductivity
th_cond_sd = 0.5          # unit: W/(mK)-1, S.D. of prior thermal conductivity
th_cond_low_bound = 0.9          # unit: W/(mK)-1, lower bound of thermal conductivity
th_cond_up_bound = 2.5          # unit: W/(mK)-1, upper bound of thermal conductivity

poro_mean = 0.3          # mean porosity
poro_sd = 0.1          # standard deviation of porosity
poro_low_bound = 0.01          # lower bound of porosity
poro_up_bound = 0.7          # upper bound of porosity