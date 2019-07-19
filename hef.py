
import numpy as np 
import pandas as pd
from src.da_methods import *
from src.pkg import *
from src.util import *
import matplotlib.pyplot as plt
import datetime
from datetime import timedelta

# Configure ES-MDA
nreaz = 100          # number of realizations
infl = 1.02          # inflation coefficient for state vector
niter = 2          # number of iterations
alpha = np.array([2,2])        # inflation coefficient for each iteration, SUM(1/alpha(i))=1 where i is the number of iteration; the coefficients could be the same or in decreasing order
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
# Configure observation 
obs_start_time = 0              # unit:s, the starting point of observation window
obs_end_time = 1203780          # unit:s, the ending point of observation window
obs_timestep = 60              # unit:s, the time interval that temperatures are collected
therm_loc = [-0.05,-0.1,-0.15] # unit:m, location of thermistor, negative means below the riverbed
obs_accuracy = 0.2              # unit: C, accuracy of temperature measurement
# Configure model domain and PFLOTRAN running environment
hz = 0.1          # unit: m, height of the 1-D column
exeprg = '/global/project/projectdirs/pflotran/pflotran-cori-new/src/pflotran/pflotran'
#exeprg = '/global/u2/k/kchen89/pflotran-cori-dbase'          # path to PFLOTRAN executable file
ncore = 100          # number of cores 

#----------------------------------------------------------
kwargs1 = {}
if 'permeability' in da_para:
    kwargs1.update({'logperm_mean':logperm_mean})
    kwargs1.update({'logperm_sd':logperm_sd})
    kwargs1.update({'perm_low_bound':perm_low_bound})
    kwargs1.update({'perm_up_bound':perm_up_bound})                 

if 'flux' in da_para:
    kwargs1.update({'flux_mean':flux_mean})
    kwargs1.update({'flux_sd':flux_sd})
    kwargs1.update({'flux_low_bound':flux_low_bound})
    kwargs1.update({'flux_up_bound':flux_up_bound}) 

if 'thermal conductivity' in da_para:
    kwargs1.update({'th_cond_mean':th_cond_mean})
    kwargs1.update({'th_cond_sd':th_cond_sd})
    kwargs1.update({'th_cond_low_bound':th_cond_low_bound})
    kwargs1.update({'th_cond_up_bound':th_cond_up_bound})
    
if 'porosity' in da_para:
    kwargs1.update({'poro_mean':poro_mean})
    kwargs1.update({'poro_sd':poro_sd})
    kwargs1.update({'poro_low_bound':poro_low_bound})
    kwargs1.update({'poro_up_bound':poro_up_bound})

th1d = TH1D(da_para,nreaz,hz,**kwargs1)
if 'permeability' in da_para:
    da_time_win = np.array([[0,obs_end_time-obs_start_time]])
else:
    da_time_win = np.array([[0,da_timestep]])
    time = da_timestep 
    while time < (obs_end_time-obs_start_time):
        da_time_win = np.append(da_time_win,[[time,time+da_timestep]],axis=0)
        time = time + da_timestep
    da_time_win = da_time_win[0:-1,:]

obs_coord = np.array(therm_loc[1:-1])-np.array(therm_loc[0])
obs_data = np.loadtxt('./observation/obs_data.dat',skiprows=1)
obs_start_idx = int(obs_start_time/obs_timestep)
obs_end_idx = int(obs_end_time/obs_timestep)
obs_data = obs_data[obs_start_idx:obs_end_idx+1,:]
obs_data[:,0] = obs_data[:,0]-obs_start_time

obs = Obs(obs_start_time,obs_end_time,obs_timestep,obs_coord,obs_accuracy,obs_data)

subprocess.call("rm -rf ./pflotran_results/*.h5",stdin=None, stdout=None,stderr=None,shell=True)
subprocess.call("rm -rf ./pflotran_results/*.chk",stdin=None, stdout=None,stderr=None,shell=True)
subprocess.call("rm -rf ./pflotran_results/*.out",stdin=None, stdout=None,stderr=None,shell=True)

subprocess.call("cp ./pflotran_inputs/1dthermal.in ./pflotran_results/",stdin=None, stdout=None,stderr=None,shell=True)

kwargs2 = {"exeprg":exeprg,"ncore":ncore,"niter":niter,"alpha":alpha,"da_time_win":da_time_win}
state_vector = Assimilator(nreaz,infl,th1d,obs,**kwargs2)