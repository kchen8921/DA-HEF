from src.pkg import *

class TH1D:
###############################################################################
#
#                Define the model parameter class TH1D
#
# Details:
#
# TH1D.state_vector:Nd X Ne, where Nd is the number of varibles and Ne is the number
#               of realizations
# TH1D.state_vector_range: Nd X 2, the first column is the lower bound for each model
#               parameter and the second column is the upper bound for each paramter
# TH1D.z: 1-D array containing the coordinates of all grid centers
#
###############################################################################
    def __init__(self,da_para,nreaz,hz,**kwargs):

        self.da_para = da_para
        self.hz = hz
        self.dz = 0.01
        self.z =(np.arange(-hz,0,self.dz)+np.arange(-hz+self.dz,self.dz,self.dz))/2

        if 'permeability' in da_para:
            logperm_mean = kwargs['logperm_mean']
            logperm_sd = kwargs['logperm_sd']
            perm_low_bound = kwargs['perm_low_bound']
            perm_up_bound = kwargs['perm_up_bound']
            init_logperm = np.random.normal(logperm_mean,logperm_sd,nreaz)
            init_perm = 10**(init_logperm)
            perm_range = np.array([perm_low_bound,perm_up_bound])
            init_perm[init_perm<perm_range[0]] = perm_range[0]
            init_perm[init_perm>perm_range[1]] = perm_range[1]
        if 'thermal conductivity' in da_para:
            th_cond_mean = kwargs['th_cond_mean']
            th_cond_sd = kwargs['th_cond_sd']
            th_cond_low_bound = kwargs['th_cond_low_bound']
            th_cond_up_bound = kwargs['th_cond_up_bound']
            init_th_cond = np.random.normal(th_cond_mean,th_cond_sd,nreaz)
            th_cond_range = np.array([th_cond_low_bound,th_cond_up_bound])
            init_th_cond[init_th_cond<th_cond_range[0]] = th_cond_range[0]
            init_th_cond[init_th_cond>th_cond_range[1]] = th_cond_range[1]
        if 'porosity' in da_para:
            poro_mean = kwargs['poro_mean']
            poro_sd = kwargs['poro_sd']
            poro_low_bound = kwargs['poro_low_bound']
            poro_up_bound = kwargs['poro_up_bound']
            init_poro = np.random.normal(poro_mean,poro_sd,nreaz)
            poro_range = np.array([poro_low_bound,poro_up_bound])
            init_poro[init_poro<poro_range[0]] = poro_range[0]
            init_poro[init_poro>poro_range[1]] = poro_range[1]
        if 'flux' in da_para:
            flux_mean = kwargs['flux_mean']
            flux_sd = kwargs['flux_sd']
            flux_low_bound = kwargs['flux_low_bound']
            flux_up_bound = kwargs['flux_up_bound']
            flux_range = np.array([flux_low_bound,flux_up_bound])
            init_flux = np.random.normal(flux_mean,flux_sd,nreaz)
            init_flux[init_flux<flux_range[0]] = flux_range[0]
            init_flux[init_flux>flux_range[1]] = flux_range[1]

        if len(da_para) == 1:
            self.state_vector = np.zeros((1,nreaz))
            self.state_vector_range = np.zeros((1,2))
            if 'permeability' in da_para:
                self.state_vector = np.array([np.log10(init_perm)])
                self.state_vector_range = np.array([np.log10(perm_range)])
            elif 'flux' in da_para:
                self.state_vector = np.array([init_flux])
                self.state_vector_range = np.array([flux_range])
            else:
                raise Exception("Please choose 'permeability' or 'flux'")
        elif len(da_para) == 2:
            self.state_vector = np.zeros((2,nreaz))
            if 'permeability' in da_para:
                self.state_vector = np.array([np.log10(init_perm)])
                self.state_vector_range = np.array([np.log10(perm_range)])
            elif 'flux' in da_para:
                self.state_vector = np.array([init_flux])
                self.state_vector_range = np.array([flux_range])
            else:
                raise Exception("Please choose 'permeability' or 'flux'")
            if 'thermal conductivity' in da_para:
                self.state_vector = np.concatenate((self.state_vector,np.array([init_th_cond])))
                self.state_vector_range = np.concatenate((self.state_vector_range,np.array([th_cond_range])))
            elif 'porosity' in da_para:
                self.state_vector = np.concatenate((self.state_vector,np.array([init_poro])))
                self.state_vector_range = np.concatenate((self.state_vector_range,np.array([poro_range])))
            else:
                raise Exception("Please choose 'thermal conductivity' or 'porosity'")
        elif len(da_para) == 3:
            self.state_vector = np.zeros((3,nreaz))
            if 'permeability' in da_para:
                self.state_vector = np.array([np.log10(init_perm)])
                self.state_vector_range = np.array([np.log10(perm_range)])
            elif 'flux' in da_para:
                self.state_vector = np.array([init_flux])
                self.state_vector_range = np.array([flux_range])
            else:
                raise Exception("Please choose 'permeability' or 'flux'")
            self.state_vector = np.concatenate((self.state_vector,np.array([init_th_cond])))
            self.state_vector = np.concatenate((self.state_vector,np.array([init_poro])))
            self.state_vector_range = np.concatenate((self.state_vector_range,np.array([th_cond_range])))
            self.state_vector_range = np.concatenate((self.state_vector_range,np.array([poro_range])))
        else:
            raise Exception("Maximum number of parameters is 3")


class Obs:
################################################################################
 #    obs: obsevation object which contains the observation time and dataself.
 #          obs.time: type: numpy.array, observation time series
 #          obs.ntime: total number of observation time points
 #          obs.value: type: numpy.array, observation data, e.g. Temperature
 #          obs.err_sd_ratio: ratio of the standard deviation of error to the observation value
 #          obs.coord: coordinates of the observation points
 #          obs.nobs: total number of observation points
################################################################################
    def __init__(self,obs_start_time,obs_end_time,obs_timestep,obs_coord,obs_accuracy,obs_data):
        self.start_time = obs_start_time
        self.end_time = obs_end_time
        self.timestep = obs_timestep
        self.time = obs_data[:,0]
        self.ntime = self.time.size
        self.value = obs_data[:,1:]
        self.coord = obs_coord
#         self.err_sd_ratio = obserr_sd_ratio
        self.accuracy = obs_accuracy
        self.nobs = obs_coord.size


def generate_dbase(nreaz,mod):
################################################################################
#
#   GenerateDbase: generate h5 file Dbase.h5 after each assimlation.
#   Dbase is a keyword in PFLOTRAN that makes the scalar value realization dependent.
#   In the TH1D model, the Dbase.h5 contains two datasets, the first is Permeability
#   and the second is ThermalConductivity. This function will be called in each iteration
#   to update the parameters.
#
#   Details:
#        nreaz: number of realizations
#        mod: model-specific object, for TH1D model, it contains the permeability
#            and thermal conductivity and associated hard limits
#
################################################################################
    filename = "./pflotran_results/Dbase.h5"
#    if os.path.isfile(filename):
#      h5file = h5py.File(filename,'r+')
#    else:
#      h5file = h5py.File(filename,'w')
    if os.path.isfile(filename):
        os.remove(filename)
        
    h5file = h5py.File(filename,'w')
    variables = []
    if 'permeability' in mod.da_para:
        variables.append("Permeability")
    elif 'flux' in mod.da_para:
        variables.append("Flux_top")
    else:
        raise Exception("Please choose 'permeability' or 'flux'")
    if 'thermal conductivity' in mod.da_para:
        variables.append('ThermalConductivity')
    if 'porosity' in mod.da_para:
        variables.append('Porosity')

    values = copy.deepcopy(mod.state_vector)
    if 'permeability' in mod.da_para:
        values[0] = 10**(values[0])
#     if 'flux' in mod.da_para:
#         values[0] = values[0]

    for i in range(len(variables)):
        if h5file.get(variables[i]):
            del h5file[variables[i]]
        h5dset = h5file.create_dataset(variables[i],data=values[i])
#    mod.state_vector[0] = np.log(mod.state_vector[0])
    h5file.close()

def generate_simulated_ensemble(nreaz,obs_coord,obs_time,z):
################################################################################
#
#   GenerateSimuEnsemble: generate the simulated data ensemble for assimlation.
#
#   Details:
#        nreaz: number of realizations
#        obs_coord: coordinates of observation points
#        obs_time: observation time
#        z: grid center of the one-dimensional model
#
################################################################################
    nobs = obs_coord.size
    obs_cell = np.zeros(nobs)
    ntime = obs_time.size
    for i in range(nobs):
        obs_cell[i] = np.argmin(np.absolute(z-obs_coord[i]))+1
    obs_cell = obs_cell.astype(int)
    simu_ensemble = np.zeros((nobs*(ntime-1),nreaz))
    for ireaz in range(nreaz):
        obs_temp = np.zeros(nobs*(ntime-1))
        j = 0
        for itime in obs_time[1:]:
            h5f = h5py.File("./pflotran_results/1dthermalR{}.h5".format(ireaz+1),'r')
            group_time = "Time:"+str(" %12.5E" % itime)+" s"
            dataset_temp = "Temperature [C]"
            obs_temp[j*nobs:(j+1)*nobs] = h5f[group_time][dataset_temp][0][0][obs_cell]
            j = j+1
        simu_ensemble[:,ireaz] = obs_temp
        #h5f.close()

    return simu_ensemble

# def run_forward_simulation_perm(nreaz,mod,obs,exeprg,ncore_per_reaz):
# ################################################################################
# #
# #   ForwardSimulationTH1D: run forward simulation of TH1D model
# #
# #   Details:
# #        nreaz: number of realizations
# #        mod: object of TH1D model
# #        obs: object of observation data
# #        exeprg: external executive program PFLOTRAN
# #        ncore: number of cores to run data assimilation on super computer
# #
# ################################################################################
#     FNULL = open(os.devnull,'w')
#     generate_dbase(nreaz,mod)
#     subprocess.call("./src/pflotran.sh {} {} {} ".format(nreaz,ncore_per_reaz,exeprg),stdin=None,stdout=FNULL,stderr=None,shell=True)
#     simu_ensemble = generate_simulated_ensemble(nreaz,obs.coord,obs.time,mod.z)
#     return simu_ensemble

# def run_forward_simulation_hy_grad(nreaz,mod,obs,exeprg,ncore_per_reaz):
# ################################################################################
# #
# #   ForwardSimulationTH1D: run forward simulation of TH1D model
# #
# #   Details:
# #        nreaz: number of realizations
# #        mod: object of TH1D model
# #        obs: object of observation data
# #        exeprg: external executive program PFLOTRAN
# #        ncore: number of cores to run data assimilation on super computer
# #
# ################################################################################
#     FNULL = open(os.devnull,'w')
#     generate_dbase(nreaz,mod)
#     make_pflotran_input(obs)
#     subprocess.call("./src/pflotran.sh {} {} {} ".format(nreaz,ncore_per_reaz,exeprg),stdin=None,stdout=FNULL,stderr=None,shell=True)
#     simu_ensemble = generate_simulated_ensemble(nreaz,obs.coord,obs.time,mod.z)
#     return simu_ensemble

def run_forward_simulation(nreaz,mod,obs,with_head,exeprg,ncore):
    if with_head:
        FNULL = open(os.devnull,'w')
        generate_dbase(nreaz,mod)
        make_pflotran_input(mod,obs,with_head,spinup=False)
        subprocess.call("./src/pflotran.sh {} {} {} ".format(nreaz,ncore,exeprg),stdin=None,stdout=FNULL,stderr=None,shell=True)

#        subprocess.call("srun -n {} {} -pflotranin ./pflotran_results/1dthermal.in -stochastic -num_realizations {} -num_groups {}".format(ncore,exeprg,nreaz,ncore))
        simu_ensemble = generate_simulated_ensemble(nreaz,obs.coord,obs.time,mod.z)
    else:
        FNULL = open(os.devnull,'w')
        generate_dbase(nreaz,mod)
        make_pflotran_input(mod,obs,with_head,spinup=False)
        subprocess.call(["./src/pflotran.sh {} {} {} ".format(nreaz,ncore,exeprg)],stdin=None,stdout=FNULL,stderr=None,shell=True)
        
#        subprocess.call("srun -n {} {} -pflotranin ./pflotran_results/1dthermal.in -stochastic -num_realizations {} -num_groups {}".format(ncore,exeprg,nreaz,ncore))
        simu_ensemble = generate_simulated_ensemble(nreaz,obs.coord,obs.time,mod.z)
#        print("obs: {} \n".format(obs.value))
#        print("mean simu_ensemble: {} \n".format(np.mean(simu_ensemble,axis=1).reshape((-1,1))))
    return simu_ensemble

def spinup(nreaz,mod,obs,with_head,exeprg,ncore):
    if with_head:
        FNULL = open(os.devnull,'w')
        generate_dbase(nreaz,mod)
        make_pflotran_input(mod,obs,with_head,spinup=True)
        subprocess.call("./src/pflotran.sh {} {} {} ".format(nreaz,ncore,exeprg),stdin=None,stdout=FNULL,stderr=None,shell=True)
#        for ireaz in range(nreaz):
#            h5f = h5py.File("./pflotran_results/1dthermalR{}.h5".format(ireaz+1),'r+')
#            group_time1 = "Time:"+str(" %12.5E" % obs.timestep)+" s"
#            group_time2 = "Time:"+str(" %12.5E" % int(259200))+" s"
#            h5f[group_time1] = h5f[group_time2]
##            del h5f[group_time2]
#            h5f.close()
    else:
        FNULL = open(os.devnull,'w')
        generate_dbase(nreaz,mod)
        make_pflotran_input(mod,obs,with_head,spinup=True)
        subprocess.call("./src/pflotran.sh {} {} {} ".format(nreaz,ncore,exeprg),stdin=None,stdout=FNULL,stderr=None,shell=True)
   
#        for ireaz in range(nreaz):
#            h5f = h5py.File("./pflotran_results/1dthermalR{}.h5".format(ireaz+1),'r+')
#            group_time1 = "Time:"+str(" %12.5E" % float(obs.timestep))+" s"
#            group_time2 = "Time:"+str(" %12.5E" % float(259200.0))+" s"
#            h5f[group_time1] = h5f[group_time2]
#            del h5f[group_time2]
#            h5f.close()

            
# def make_pflotran_input(obs):
#     obs_win = np.loadtxt('./observation/da_time_win.dat')
#     with open('./pflotran_results/1dthermal.in','r+') as f:
#         pflotranin = f.readlines()
        
#         if (obs.time[0]-obs_win[0]) > 1e-3:
#             restart_lindex = [i for i,s in enumerate(pflotranin) if "FILENAME 1dthermal" in s][0]
#             pflotranin[restart_lindex-1] = "  RESTART"+"\n"
#             pflotranin[restart_lindex] = "    FILENAME 1dthermal-restart.chk "+"\n"
#             pflotranin[restart_lindex+1] = "    REALIZATION_DEPENDENT"+"\n"
#         final_time_lindex = [i for i,s in enumerate(pflotranin) if "FINAL_TIME" in s][0]
#         pflotranin[final_time_lindex] = "  FINAL_TIME {} sec".format(obs.time[-1])+"\n"
#     f.close()
#     os.remove('./pflotran_results/1dthermal.in')
#     with open('./pflotran_results/1dthermal.in','w') as new_f:
#         new_f.writelines(pflotranin)
#     new_f.close()
#     return

def initialize_restart_file(nreaz,obs_time):
    for ireaz in range(nreaz):
        h5f = h5py.File("./pflotran_results/1dthermalR{}.h5".format(ireaz+1),'r+')
        for itime in obs_time[1:]:
            group_time = "Time:"+str(" %12.5E" % itime)+" s"
            del h5f[group_time]
        h5f.close()
        
def copy_chk():
    subprocess.call("cp ./pflotran_results/*restart.chk ./pflotran_results/chk_temp",stdin=None, stdout=None,stderr=None,shell=True)
    
def load_chk():
    subprocess.call("cp ./pflotran_results/chk_temp/*restart.chk ./pflotran_results/",stdin=None, stdout=None,stderr=None,shell=True)
    
def make_pflotran_input(mod,obs,with_head,spinup):
#    obs_win = np.loadtxt('./observation/da_time_win.dat')
    with open('./pflotran_results/1dthermal.in','r+') as f:
        pflotranin = f.readlines()
        if with_head:
            for i,s in enumerate(pflotranin):
                if "PERM_ISO" in s:
                    pflotranin[i] = "    PERM_ISO DBASE_VALUE Permeability" + "\n"
#                 if 'TEMPERATURE FILE ../pflotran_inputs/temp_top.dat' in s:
#                     pflotranin[i-3] = "  DATUM FILE ../pflotran_inputs/head_top.dat" + "\n"
#                 if 'TEMPERATURE FILE ../pflotran_inputs/temp_bottom.dat' in s:
#                     pflotranin[i-2] = "  DATUM FILE ../pflotran_inputs/head_bottom.dat" + "\n"                             
                if 'NXYZ' in s:
                    pflotranin[i] = "  NXYZ 1 1 {}".format(int(mod.hz*100)) + "\n"
                    pflotranin[i+2] = "    0.d0 0.d0 {}".format(-mod.hz) + "d0" + "\n"
                if 'REGION all' in s and 'COORDINATES' in pflotranin[i+1]:
                    pflotranin[i+2] = "    0.d0 0.d0 {}".format(-mod.hz) + "d0" + "\n"
                if 'REGION bottom' in s and "FACE" in pflotranin[i+1]:
                    pflotranin[i+3] = "    0.d0 0.d0 {}".format(-mod.hz) + "d0" + "\n"
                    pflotranin[i+4] = "    1.d0 1.d0 {}".format(-mod.hz) + "d0" + "\n"
                if 'SNAPSHOT_FILE' in s:
                    if spinup:
                        pflotranin[i+1] = "   PERIODIC TIME {}".format(259200) + " s" +"\n"
                    else:
                        pflotranin[i+1] = "   PERIODIC TIME {}".format(obs.timestep) + " s" +"\n"
                if 'FLOW_CONDITION flow_top' in s and "TYPE" in pflotranin[i+1]:
                    pflotranin[i+5] = "  DATUM FILE ../pflotran_inputs/head_top.dat" + "\n"
                if 'FLOW_CONDITION flow_bottom' in s and "TYPE" in pflotranin[i+1]:
                    pflotranin[i+5] = "  DATUM FILE ../pflotran_inputs/head_bottom.dat" + "\n"
                if 'FLOW_CONDITION initial' in s and "TYPE" in pflotranin[i+1]:
                    pflotranin[i+6] = "  TEMPERATURE " + str(np.mean(obs.value[0,:])) + "d0" + "\n"
                if 'THERMAL_CONDUCTIVITY_WET' in s:
                    if 'thermal conductivity' in mod.da_para:
                        pflotranin[i] = "  THERMAL_CONDUCTIVITY_WET DBASE_VALUE ThermalConductivity"
                if 'POROSITY' in s:
                    if 'porosity' in mod.da_para:
                        pflotranin[i] = "  POROSITY DBASE_VALUE Porosity"
                if "FILENAME 1dthermal" in s:
                    if not spinup:
                        pflotranin[i-1] = "  RESTART"+"\n"
                        pflotranin[i] = "    FILENAME 1dthermal-restart.chk "+"0.d0 \n"
                        pflotranin[i+1] = "    REALIZATION_DEPENDENT"+"\n"    
                        pflotranin[i+2] = "    RESET_TO_TIME_ZERO" +"\n"                     
                if 'FINAL_TIME' in s:
                    if spinup:
                        pflotranin[i] = "  FINAL_TIME " + str(86400) + "  sec" + "\n"
                    else:       
                        pflotranin[i] = "  FINAL_TIME " + str(obs.time[-1]) + "  sec" + "\n"
                    
        else:
            for i,s in enumerate(pflotranin):
                if 'NXYZ' in s:
                    pflotranin[i] = "  NXYZ 1 1 {}".format(int(mod.hz*100)) + "\n"
                    pflotranin[i+2] = "    0.d0 0.d0 {}".format(-mod.hz) + "d0" + "\n"
                if 'REGION all' in s and 'COORDINATES' in pflotranin[i+1]:
                    pflotranin[i+2] = "    0.d0 0.d0 {}".format(-mod.hz) + "d0" + "\n"
                if 'REGION bottom' in s and "FACE" in pflotranin[i+1]:
                    pflotranin[i+3] = "    0.d0 0.d0 {}".format(-mod.hz) + "d0" + "\n"
                    pflotranin[i+4] = "    1.d0 1.d0 {}".format(-mod.hz) + "d0" + "\n"
                if 'SNAPSHOT_FILE' in s:
                    if spinup:
                        pflotranin[i+1] = "   PERIODIC TIME {}".format(259200) + " s" +"\n"
                    else:
                        pflotranin[i+1] = "    PERIODIC TIME {}".format(obs.timestep) + " s" +"\n"
                if "FLOW_CONDITION flow_top" in s and "TYPE" in pflotranin[i+1]:
                    pflotranin[i+2] = "    FLUX NEUMANN" + "\n"
                    pflotranin[i+5] = "\n"
                    pflotranin[i+6] = "  FLUX DBASE_VALUE Flux_top m/day" +"\n"
                if 'FLOW_CONDITION initial' in s and "TYPE" in pflotranin[i+1]:
                    pflotranin[i+6] = "  TEMPERATURE " + str(np.mean(obs.value[0,:])) + "d0" + "\n"
                if 'THERMAL_CONDUCTIVITY_WET' in s:
                    if 'thermal conductivity' in mod.da_para:
                        pflotranin[i] = "  THERMAL_CONDUCTIVITY_WET DBASE_VALUE ThermalConductivity" + "\n"
                if 'POROSITY' in s:
                    if 'porosity' in mod.da_para:
                        pflotranin[i] = "  POROSITY DBASE_VALUE Porosity" + "\n"
                if "FINAL_TIME" in s:
                    pflotranin[i] = "  FINAL_TIME {} sec".format(obs.time[-1])+"\n"

                if "FILENAME 1dthermal" in s:
                    if not spinup:
#                        print(pflotranin[i+2])
                        pflotranin[i-1] = "  RESTART"+"\n"
                        pflotranin[i] = "    FILENAME 1dthermal-restart.chk "+" \n"
                        pflotranin[i+1] = "    REALIZATION_DEPENDENT"+"\n"  
                        if obs.time[0] < 2e4:   
                          pflotranin[i+2] = "    RESET_TO_TIME_ZERO" +"\n"   
                        else: 
                          pflotranin[i+2] = "#    RESET_TO_TIME_ZERO" +"\n"
                if "FINAL_TIME" in s:
                    if spinup:
                        pflotranin[i] = "  FINAL_TIME {} sec".format(86400)+"\n"
                    else:
                        pflotranin[i] = "  FINAL_TIME {} sec".format(obs.time[-1])+"\n"                 
    f.close()

    os.remove('./pflotran_results/1dthermal.in')
    with open('./pflotran_results/1dthermal.in','w') as new_f:
        new_f.writelines(pflotranin)
    new_f.close()
    return    

def draw_plot(data, edge_color, fill_color1,fill_color2):
    bp = plt.boxplot(data, patch_artist=True, showfliers=False)
    colors = [fill_color1,fill_color2,fill_color2,fill_color2,fill_color2]
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)

    for patch,color in zip(bp['boxes'],colors):
        patch.set(facecolor=color)
        
def plot_para_without_head(da_para,da_timestep,obs_timestep,nreaz):
    state_vector = np.loadtxt('./pflotran_results/state_vector_out_without_head.txt')
    obs_data = np.loadtxt('./pflotran_results/obs_out_without_head.txt')

    
    sec_to_day = 3600*24
    obs_data[:,0] = obs_data[:,0]/sec_to_day

    with open('./pflotran_inputs/1dthermal.in','r+') as f:
        pflotranin = f.readlines()
        for i,s in enumerate(pflotranin):
            if 'PERM_ISO' in s:
                perm = float(s.split()[1])
    f.close()


        
    fig = plt.figure(num=1,dpi=300)
    fig.subplots_adjust(hspace=0.3,wspace=0.3)
    plt.tight_layout()

    state_vector = state_vector[len(da_para):,:]
    for i in range(len(da_para)):
        plt.subplot(len(da_para),1,i+1)
        if i == 0:
            m = int(da_timestep/obs_timestep)

            n = int(state_vector.shape[0]/len(da_para))
            obs_time = np.arange(0,n*da_timestep,da_timestep)/sec_to_day
            line1, = plt.plot(obs_time,state_vector[:,0],'chocolate',linewidth=0.3)
            for j in range(nreaz-1):
                plt.plot(obs_time,state_vector[:,j+1],'chocolate',linewidth=0.3)
            line2, = plt.plot(obs_time,np.mean(state_vector,axis=1),'b',linewidth=1)
            line3, = plt.plot(obs_time,np.percentile(state_vector,95,axis=1),'b--',linewidth=0.5)
            line4, = plt.plot(obs_time,np.percentile(state_vector,5,axis=1),'b--',linewidth=0.5)
            
            if len(da_para) == 1:
                plt.xlabel('Time (day)')
            plt.ylabel('q (m/d)')
            ymax = np.ceil(max(np.mean(hef_vector,axis=1)))
            ymin = np.floor(min(np.mean(hef_vector,axis=1)))
            nticks = 4
            plt.yticks(np.arange(ymin, ymax+0.1, (ymax-ymin)/nticks))
            plt.legend((line1,line2,line3,line4),('Posterior','Mean','95%','5%'),frameon=False,ncol=4)    
        else:
            plt.subplot(len(da_para),1,i+1)
            line1, = plt.plot(obs_time,state_vector[i::len(da_para),0],'chocolate',linewidth=0.3)
            for j in range(nreaz-1):
                plt.plot(obs_time,state_vector[i::len(da_para),j+1],'chocolate',linewidth=0.3)
            line2, = plt.plot(obs_time,np.mean(state_vector[i::len(da_para),:],axis=1),'b',linewidth=1)
            line3, = plt.plot(obs_time,np.percentile(state_vector[i::len(da_para),:],95,axis=1),'b--',linewidth=0.5)
            line4, = plt.plot(obs_time,np.percentile(state_vector[i::len(da_para),:],5,axis=1),'b--',linewidth=0.5)

            if da_para[i] == 'thermal conductivity':
                plt.ylabel('$\lambda$ (W/m$\cdot$K)')
            else:
                plt.ylabel('$\phi$')
            ymax = np.ceil(max(np.mean(state_vector[i::len(da_para),:],axis=1)))
            ymin = np.floor(min(np.mean(state_vector[i::len(da_para),:],axis=1)))
            nticks = 4
            plt.yticks(np.arange(ymin, ymax+0.1, (ymax-ymin)/nticks))

            if i == len(da_para)-1:
                plt.xlabel('Time (day)')
                
def plot_para_with_head(da_para):
    state_vector = np.loadtxt('./pflotran_results/state_vector_out_with_head.txt')
    obs_data = np.loadtxt('./pflotran_results/obs_out_with_head.txt')
    
    sec_to_day = 3600*24
    obs_data[:,0] = obs_data[:,0]/sec_to_day

    fig = plt.figure(num=1,dpi=300)
    fig.subplots_adjust(wspace=0.3)
    plt.tight_layout()

    for i in range(len(da_para)):
        plt.subplot(1,len(da_para),i+1)
        perm = np.array([state_vector[i,:],state_vector[i+len(da_para),:]])
        draw_plot(np.transpose(perm),'red','cyan','tan')
        if da_para[i] == 'permeability':
            plt.ylabel('log$_{10}(k$) (m$^{2}$)')
        elif da_para[i] == 'thermal conductivity':
            plt.ylabel('Thermal Cond. (W/m$\cdot$K)')
        else:
            plt.ylabel('Porosity')
            
        plt.title('{}'.format(da_para[i]))
        plt.xticks([1,2],['Prior','Posterior'])
        
def plot_temp_without_head(nreaz,therm_loc,da_para,obs_timestep,init_datetime):
    state_vector = np.loadtxt('./pflotran_results/state_vector_out_without_head.txt')
    simu_ensemble = np.loadtxt('./pflotran_results/simu_ensemble_out_without_head.txt')
    obs_data = np.loadtxt('./pflotran_results/obs_out_without_head.txt')

    sec_to_day = 3600*24
    obs_data[:,0] = obs_data[:,0]/sec_to_day
    nobs = len(therm_loc)-2
    init_datetime = datetime.datetime.strptime(init_datetime,"%Y-%m-%d %H:%M:%S")

    obs_time1 = []
    for i in range(obs_data.shape[0]):
        obs_time1.append(init_datetime+timedelta(seconds=obs_data[i,0]))

    fig = plt.figure(num=1,dpi=150)
    for i in range(nobs):
        plt.subplot(nobs,1,i+1)
        line1, = plt.plot(obs_time1,obs_data[:,i+1],'k',linewidth=1)
        line2, = plt.plot(obs_time1,simu_ensemble[i::nobs,0],'chocolate',linewidth=0.3)
        for j in range(nreaz-1):
            plt.plot(obs_time1,simu_ensemble[i::nobs,j+1],'chocolate',linewidth=0.3)
        line3, = plt.plot(obs_time1,np.mean(simu_ensemble[i::nobs,:],axis=1),'b',linewidth=1)


        plt.xlabel('Time (day)')
        plt.ylabel('Temperature ($^\circ$C)')
        plt.legend((line1,line2,line3),('Observation','Posterior','Mean'),frameon=False)
        plt.title('Posterior: Obs. point {} m'.format(therm_loc[i+1]))  
        ax = plt.gca()
        dayFmt = mdates.DateFormatter("%m/%d/%y")
        days = mdates.DayLocator(interval=3)    
        ax.xaxis.set_major_locator(days)
        ax.xaxis.set_major_formatter(dayFmt)
        ax.xaxis.set_minor_locator(mdates.DayLocator())
        
def plot_temp_with_head(obs_coord,nreaz):
    state_vector = np.loadtxt('./pflotran_results/state_vector_out_with_head.txt')
    simu_ensemble = np.loadtxt('./pflotran_results/simu_ensemble_out_with_head.txt')
    simu_ensemble_prior = np.loadtxt('./pflotran_results/simu_ensemble_out_with_head_prior.txt')
    obs_data = np.loadtxt('./pflotran_results/obs_out_with_head.txt')

    sec_to_day = 3600*24
    obs_data[:,0] = obs_data[:,0]/sec_to_day
    fig = plt.figure(num=1,dpi=150)
    nobs = len(obs_coord)
    
    for i in range(nobs):
        plt.subplot(nobs,2,i*2+1)
        line1, = plt.plot(obs_data[:,0],simu_ensemble_prior[i::nobs,0],'chocolate',linewidth=0.3)
        for j in range(nreaz-1):
            plt.plot(obs_data[:,0],simu_ensemble_prior[i::nobs,j+1],'chocolate',linewidth=0.3)
        line2, = plt.plot(obs_data[:,0],np.mean(simu_ensemble_prior[i::nobs,:],axis=1),'b',linewidth=1)
        line3, = plt.plot(obs_data[:,0],obs_data[:,i+1],'k',linewidth=1)
        plt.xlabel('Time (day)')
        plt.ylabel('Temperature ($^\circ$C)')
        plt.legend((line1,line2,line3),('Posterior','Mean','Observation'),frameon=False)
        plt.title('Prior: Obs. point {}'.format(obs_coord[i]))

        plt.subplot(nobs,2,i*2+2)
        line1, = plt.plot(obs_data[:,0],simu_ensemble[i::nobs,0],'chocolate',linewidth=0.3)
        for j in range(nreaz-1):
            plt.plot(obs_data[:,0],simu_ensemble[i::nobs,j+1],'chocolate',linewidth=0.3)
        line2, = plt.plot(obs_data[:,0],np.mean(simu_ensemble[i::nobs,:],axis=1),'b',linewidth=1)
        line3, = plt.plot(obs_data[:,0],obs_data[:,i+1],'k',linewidth=1)

        plt.xlabel('Time (day)')
#         plt.ylabel('Temperature ($^\circ$C)')
        plt.legend((line1,line2,line3),('Posterior','Mean','Observation'),frameon=False)
        plt.title('Posterior: Obs. point {}'.format(obs_coord[i]))
        
def plot_hef_without_head(da_para,da_timestep,obs_timestep,nreaz,init_datetime): 
    state_vector = np.loadtxt('./pflotran_results/state_vector_out_without_head.txt')
    simu_ensemble = np.loadtxt('./pflotran_results/simu_ensemble_out_without_head.txt')
    obs_data = np.loadtxt('./pflotran_results/obs_out_without_head.txt')

    sec_to_day = 3600*24
    with open('./pflotran_inputs/1dthermal.in','r+') as f:
        pflotranin = f.readlines()
        for i,s in enumerate(pflotranin):
            if 'PERM_ISO' in s:
                perm = float(s.split()[1])
    f.close()

    init_datetime = datetime.datetime.strptime(init_datetime,"%Y-%m-%d %H:%M:%S")
    state_vector = state_vector[len(da_para):,:]
    m = int(da_timestep/obs_timestep)
    hef_vector = state_vector[::len(da_para),:]

    fig = plt.figure(num=1,dpi=150)

    n = int(state_vector.shape[0]/len(da_para))
    obs_time = np.arange(0,n*da_timestep,da_timestep)
    obs_time1 = []
    for i in range(obs_time.shape[0]):
        obs_time1.append(init_datetime+timedelta(seconds=obs_time[i]))

    line1, = plt.plot(obs_time1,hef_vector[:,0],'chocolate',linewidth=0.3)
    for j in range(nreaz-1):
        plt.plot(obs_time1,hef_vector[:,j+1],'chocolate',linewidth=0.3)
    line2, = plt.plot(obs_time1,np.mean(hef_vector,axis=1),'b',linewidth=1)
    line3, = plt.plot(obs_time1,np.percentile(hef_vector,95,axis=1),'b--',linewidth=0.5)
    line4, = plt.plot(obs_time1,np.percentile(hef_vector,5,axis=1),'b--',linewidth=0.5)
    line5, = plt.plot(obs_time1,np.zeros(len(obs_time)),'k--',linewidth=0.3)
    plt.xlabel('Time (day)')
    plt.ylabel('HEF (m/d)')
    plt.legend((line1,line2,line3,line4),('Posterior','Mean','95%','5%'),frameon=False) 
    ax = plt.gca()
    dayFmt = mdates.DateFormatter("%m/%d/%y")
    days = mdates.DayLocator(interval=3)    
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(dayFmt)
    ax.xaxis.set_minor_locator(mdates.DayLocator())
    plt.savefig('./pflotran_results/hef_without_head.png',dpi=150)
    
def plot_hef_with_head(da_para,hz,obs_start_time,obs_end_time,nreaz):
    state_vector = np.loadtxt('./pflotran_results/state_vector_out_with_head.txt')
    simu_ensemble = np.loadtxt('./pflotran_results/simu_ensemble_out_with_head.txt')
    obs_data = np.loadtxt('./pflotran_results/obs_out_with_head.txt')
    head_top = np.loadtxt('./pflotran_inputs/head_top.dat')
    head_bottom = np.loadtxt('./pflotran_inputs/head_bottom.dat')
    obs_start_idx = np.argmin(np.absolute(obs_start_time-head_top[:,0]))
    obs_end_idx = np.argmin(np.absolute(obs_end_time-head_top[:,0]))
    head_top = head_top[obs_start_idx:obs_end_idx,:]
    head_bottom = head_bottom[obs_start_idx:obs_end_idx,:]

    sec_to_day = 3600*24
    perm = 10**(state_vector[len(da_para),:])
    temp = obs_data[:,1]+273.15
    hy_grad = (head_top[:,3]-head_bottom[:,3])/hz
    viscosity = 1e-6*(280.68*(temp/300)**(-1.9)+511.45*(temp/300)**(-7.7)+61.131*(temp/300)**(-19.6)+0.45903*(temp/300)**(-40))
    hef_vector = sec_to_day*1000*9.8*np.matmul(np.multiply(np.transpose(hy_grad),1/viscosity).reshape(-1,1),perm.reshape(1,-1))

    fig = plt.figure(num=1,dpi=150)

    n = int(state_vector.shape[0]/len(da_para))
    obs_time = obs_data[:,0]/sec_to_day
    line1, = plt.plot(obs_time,hef_vector[:,0],'chocolate',linewidth=0.3)
    for j in range(nreaz-1):
        plt.plot(obs_time,hef_vector[:,j+1],'chocolate',linewidth=0.3)
    line2, = plt.plot(obs_time,np.mean(hef_vector,axis=1),'b',linewidth=1)
    line3, = plt.plot(obs_time,np.percentile(hef_vector,95,axis=1),'b--',linewidth=0.5)
    line4, = plt.plot(obs_time,np.percentile(hef_vector,5,axis=1),'b--',linewidth=0.5)
    line5, = plt.plot(obs_time,np.zeros(len(obs_time)),'k--',linewidth=0.3)

    plt.xlabel('Time (day)')
    plt.ylabel('HEF (m/d)')
    plt.legend((line1,line2,line3,line4),('Posterior','Mean','95%','5%'),frameon=False)         
    plt.savefig('./pflotran_results/hef_with_head.png',dpi=150)

def plot_temp_data(therm_loc,obs_start_time,obs_end_time,init_datetime):
    temp_top = np.loadtxt('./pflotran_inputs/temp_top.dat',skiprows=1)
    temp_bot = np.loadtxt('./pflotran_inputs/temp_bottom.dat',skiprows=1)
    temp_obs = np.loadtxt('./observation/obs_data.dat',skiprows=1)
    init_datetime = datetime.datetime.strptime(init_datetime,"%Y-%m-%d %H:%M:%S")

    sec_to_day = 3600*24
    start_idx = np.argmin(np.abs(temp_obs[:,0]-obs_start_time))
    end_idx = np.argmin(np.abs(temp_obs[:,0]-obs_end_time))
    time = temp_obs[start_idx:end_idx,0]/sec_to_day
    time1 = []
    for i in range(time.shape[0]):
        time1.append(init_datetime+timedelta(days=time[i]))

    temp_top = temp_top[start_idx:end_idx,:]
    temp_bot = temp_bot[start_idx:end_idx,:]
    temp_obs = temp_obs[start_idx:end_idx,:]

    fig = plt.figure(num=1,dpi=150)

    line1, = plt.plot(time1,temp_top[:,1],'b',linewidth=1)
    line2, = plt.plot(time1,temp_bot[:,1],'r',linewidth=1)

    lines = [line1,line2]
    legends = [str(therm_loc[0])+' m(top)',str(therm_loc[-1])+' m(bottom)']
    colors = ['g','c','k','y']
    for i in range(len(therm_loc)-2):
        line, = plt.plot(time1,temp_obs[:,i+1],linewidth=1,color=colors[i])
        lines.append(line)
        legends.append(str(therm_loc[i+1])+' m')

    plt.xlabel('Time (day)')
    plt.ylabel('Temperature ($^\circ$C)')
    plt.legend(lines,legends,frameon=False)
    plt.title('Temperature data')
    ax = plt.gca()
    dayFmt = mdates.DateFormatter("%m/%d/%y")
    days = mdates.DayLocator(interval=3)    
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(dayFmt)
    ax.xaxis.set_minor_locator(mdates.DayLocator())
    plt.savefig('./pflotran_results/temperature_data.png',dpi=150)