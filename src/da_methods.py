from src.pkg import *
from src.util import *

# def Assimilator(nreaz,infl,mod,obs,continuous,**kwargs):
#     state_vector_out = mod.state_vector

# ################################################################################
# #                  Ensemble Smoother-Multiple Data Assimilation
# #
# # reference: Emerick, A., and Reynolds, A.,2013, Ensemble smoother with Multiple
# # data assimilation, Comput. Geosci.
# #
# # Paramters:
# #    nreaz: number of realizations
# #    infl: inflation coefficients
# #    mod: TH1D object that contains all the model related parameters
# #    obs: observation object that contains all the observation data
# #    **kwargs:
# #          "exeprg": path to the the executable program
# #          "ncore_per_reaz": number of cores needed for running single realization
# #          "niter": number of iterations for ensember smoother
# #          "alpha": inflation coefficients for observation error covariance. The requirements
# #                   for alpha is SUM(1/alpha(i)) = 1 where i is the number of iteration.
# #          "da_time_win": time interval in which data are used for assimilation
# #
# # Return:
# #        state_vector_out: updated state vector after every assimilation is appended to
# #                          this vector
# ###############################################################################

#     if continuous == False:
#         obs_flow_dir = np.loadtxt('./observation/obs_flow_dir.dat')

#     if len(kwargs)!=5:
#         raise Exception("Please provide five arguements following the order: exeprg, ncore_per_reaz, niter, alpha, ndata")
#     elif len(kwargs['alpha'])!=kwargs['niter']:
#         raise Exception("length of alpha should be equal to niter")
#     else:
#         exeprg = kwargs['exeprg']
#         ncore_per_reaz = kwargs['ncore_per_reaz']
#         niter = kwargs['niter']
#         alpha = kwargs['alpha']
#         da_time_win = kwargs['da_time_win']

#     f = open("./results.txt",'w')
#     ibatch = 0 # one batch corresponds to one data assimilation time interval
#     while ibatch<len(da_time_win):
#         if continuous == False and ibatch>0:
#             if obs_flow_dir[ibatch][1]>0 and obs_flow_dir[ibatch-1][1]<0:
#                 mod.state_vector[:,:] = np.random.normal(0.02,0.01,nreaz)
#             elif obs_flow_dir[ibatch][1]<0 and obs_flow_dir[ibatch-1][1]>0:
#                 mod.state_vector[:,:] = np.random.normal(-0.02,0.01,nreaz)

#         print("time win is {} \n".format(da_time_win[ibatch]))
#         f.write("time window is {} \n".format(da_time_win[ibatch]))
#         idx_time_start = np.argmin(abs(obs.time-da_time_win[ibatch][0]))
#         idx_time_end = np.argmin(abs(obs.time-da_time_win[ibatch][1]))
#         obs_temp = copy.deepcopy(obs)
#         obs_temp.time = obs_temp.time[idx_time_start:idx_time_end+1]
#         obs_temp.value = obs_temp.value[idx_time_start:idx_time_end+1]
#         obs_temp.ntime = len(obs_temp.time)
#         for i in range(niter):
#             if 'permeability' in mod.da_para:
#             # prediction
#                 simu_ensemble = run_forward_simulation_perm(nreaz,mod,obs_temp,exeprg,ncore_per_reaz)
#             else:
#                 simu_ensemble = run_forward_simulation_hy_grad(nreaz,mod,obs_temp,exeprg,ncore_per_reaz)
                
#             # analysis
#             obs_sd = math.sqrt(alpha[i])*obs_temp.err_sd_ratio*obs_temp.value
#             obs_ensemble = np.repeat(obs_temp.value.flatten('C').reshape(obs_temp.ntime*obs_temp.nobs,1),nreaz,1)+np.dot(np.diag(obs_sd.flatten('C')),np.random.normal(0,1,nreaz*obs_temp.nobs*obs_temp.ntime).reshape(obs_temp.nobs*obs_temp.ntime,nreaz))
#             state_vector = mod.state_vector[:,:]
#             cov_state_simu = np.cov(state_vector,simu_ensemble)[0:len(state_vector),len(state_vector):]
#             cov_simu = np.cov(simu_ensemble)
#             if obs.nobs*obs.ntime == 1:
#                 inv_cov_simu_obserr = np.array(1/(cov_simu+np.square(np.diag(obs_sd))))
#             else:
#                 inv_cov_simu_obserr = la.inv(cov_simu+np.square(np.diag(obs_sd.flatten('C'))))
#             kalman_gain = np.dot(cov_state_simu,inv_cov_simu_obserr)
#             state_vector = state_vector+np.dot(kalman_gain,obs_ensemble-simu_ensemble)
#             print("mean is {} ".format(np.mean(state_vector[0])))
#             print("std is {} ".format(np.std(state_vector[0])))
#             for j in range(state_vector.shape[0]):
#                 state_vector[j,:][state_vector[j,:]<mod.state_vector_range[j,0]] = mod.state_vector_range[j,0]
#                 state_vector[j,:][state_vector[j,:]>mod.state_vector_range[j,1]] = mod.state_vector_range[j,1]
#             mod.state_vector[:,:] = state_vector
#         state_vector_out = np.concatenate((state_vector_out,state_vector),axis=0)
#         np.savetxt('./results/state_vector_out.txt',state_vector_out)
#         mod.state_vector = inflate_state_vector(state_vector,infl)
#         ibatch = ibatch+1
#     return state_vector_out

def Assimilator(nreaz,infl,mod,obs,**kwargs):

################################################################################
#                  Ensemble Smoother-Multiple Data Assimilation
#
# reference: Emerick, A., and Reynolds, A.,2013, Ensemble smoother with Multiple
# data assimilation, Comput. Geosci.
#
# Paramters:
#    nreaz: number of realizations
#    infl: inflation coefficients
#    mod: TH1D object that contains all the model related parameters
#    obs: observation object that contains all the observation data
#    **kwargs:
#          "exeprg": path to the the executable program
#          "ncore_per_reaz": number of cores needed for running single realization
#          "niter": number of iterations for ensember smoother
#          "alpha": inflation coefficients for observation error covariance. The requirements
#                   for alpha is SUM(1/alpha(i)) = 1 where i is the number of iteration.
#          "da_time_win": time interval in which data are used for assimilation
#
# Return:
#        state_vector_out: updated state vector after every assimilation is appended to
#                          this vector
###############################################################################

#     if continuous == False:
#         obs_flow_dir = np.loadtxt('./observation/obs_flow_dir.dat')

    if 'permeability' in mod.da_para:
        with_head = True
    else:
        with_head = False
            
    if len(kwargs)!=5:
        raise Exception("Please provide five arguements following the order: exeprg, ncore_per_reaz, niter, alpha, ndata")
    elif len(kwargs['alpha'])!=kwargs['niter']:
        raise Exception("length of alpha should be equal to niter")
    else:
        exeprg = kwargs['exeprg']
        ncore = kwargs['ncore']
        niter = kwargs['niter']
        alpha = kwargs['alpha']
        da_time_win = kwargs['da_time_win']

    spinup(nreaz,mod,obs,with_head,exeprg,ncore)

#     f = open("./results.txt",'w')
    state_vector_out = copy.deepcopy(mod.state_vector)
    ibatch = 0 # one batch corresponds to one data assimilation time interval
    while ibatch<len(da_time_win):
#         if continuous == False and ibatch>0:
#             if obs_flow_dir[ibatch][1]>0 and obs_flow_dir[ibatch-1][1]<0:
#                 mod.state_vector[:,:] = np.random.normal(0.02,0.01,nreaz)
#             elif obs_flow_dir[ibatch][1]<0 and obs_flow_dir[ibatch-1][1]>0:
#                 mod.state_vector[:,:] = np.random.normal(-0.02,0.01,nreaz)

        print("Data assimilation time window: {} ".format(da_time_win[ibatch]))
#         f.write("time window is {} \n".format(da_time_win[ibatch]))
        idx_time_start = np.argmin(abs(obs.time-da_time_win[ibatch][0]))
        idx_time_end = np.argmin(abs(obs.time-da_time_win[ibatch][1]))
        obs_temp = copy.deepcopy(obs)
        obs_temp.time = obs_temp.time[idx_time_start:idx_time_end+1]
        obs_temp.value = obs_temp.value[idx_time_start:idx_time_end+1]
        obs_temp.ntime = len(obs_temp.time)
        for i in range(niter):
            if ibatch > 0 and i > 0:
                load_chk()
                initialize_restart_file(nreaz,obs_temp.time)
#                time.sleep(300)
                
            if ibatch > 0 and i == 0:
                copy_chk()            
                                
            simu_ensemble = run_forward_simulation(nreaz,mod,obs_temp,with_head,exeprg,ncore)
#            print("simu_ensemble: ",simu_ensemble[1,:])
            if with_head and i == 0:
                np.savetxt('./pflotran_results/simu_ensemble_out_with_head_prior.txt',simu_ensemble[0:obs_temp.ntime-1,:])
#             if 'permeability' in mod.da_para:
#             # prediction
#                 simu_ensemble = run_forward_simulation_perm(nreaz,mod,obs_temp,exeprg,ncore_per_reaz)
#             else:
#                 simu_ensemble = run_forward_simulation_hy_grad(nreaz,mod,obs_temp,exeprg,ncore_per_reaz)
                
            # analysis
#             obs_sd = math.sqrt(alpha[i])*obs_temp.err_sd_ratio*obs_temp.value
#             obs_ensemble = np.repeat(obs_temp.value.flatten('C').reshape(obs_temp.ntime*obs_temp.nobs,1),nreaz,1)+np.dot(np.diag(obs_sd.flatten('C')),np.random.normal(0,1,nreaz*obs_temp.nobs*obs_temp.ntime).reshape(obs_temp.nobs*obs_temp.ntime,nreaz))
            obs_err_sd = obs_temp.accuracy/3*np.ones(obs_temp.value[1:].shape)
            obs_ensemble = np.repeat(obs_temp.value[1:].flatten('C').reshape((obs_temp.ntime-1)*obs_temp.nobs,1),nreaz,1)+np.sqrt(alpha[i])*np.dot(np.diag(obs_err_sd.flatten('C')),np.random.normal(0,1,nreaz*obs_temp.nobs*(obs_temp.ntime-1)).reshape(obs_temp.nobs*(obs_temp.ntime-1),nreaz))
            state_vector = mod.state_vector[:,:]
            cov_state_simu = np.cov(state_vector,simu_ensemble)[0:len(state_vector),len(state_vector):]
            cov_simu = np.cov(simu_ensemble)

            if obs.nobs*obs.ntime == 1:
                inv_cov_simu_obserr = np.array(1/(cov_simu+np.square(np.diag(obs_err_sd))))
            else:
                inv_cov_simu_obserr = la.inv(cov_simu+np.square(np.diag(obs_err_sd.flatten('C'))))
            kalman_gain = np.dot(cov_state_simu,inv_cov_simu_obserr)
            state_vector = state_vector+np.dot(kalman_gain,obs_ensemble-simu_ensemble)
            
            if with_head:
                print("Iteration {}:".format(i))
                print("  Mean of log(permeability) is {} ".format(np.mean(state_vector[0])))
                print("  STD of log(permeability) is {} \n".format(np.std(state_vector[0])))
            else:
                if i == niter-1:
                    print("  Mean of flux after {} iterations is {} ".format(niter,np.mean(state_vector[0])))
                    print("  STD of flux after {} iterations is {} \n ".format(niter,np.std(state_vector[0])))                    
            for j in range(state_vector.shape[0]):
                state_vector[j,:][state_vector[j,:]<mod.state_vector_range[j,0]] = mod.state_vector_range[j,0]
                state_vector[j,:][state_vector[j,:]>mod.state_vector_range[j,1]] = mod.state_vector_range[j,1]
            mod.state_vector[:,:] = state_vector
#             if with_head:
#                 np.savetxt('./pflotran_results/state_vector_iter{}.txt'.format(i+1),state_vector)


            
        n = (obs_temp.ntime-1)*obs_temp.nobs
        if ibatch == 0:            
            simu_ensemble_out = simu_ensemble[0:obs_temp.ntime-1,:]
            obs_out = np.concatenate((obs_temp.time[0:obs_temp.ntime-1].reshape(-1,1),obs_temp.value[0:obs_temp.ntime-1,:]),axis=1)
        else:
            simu_ensemble_out = np.concatenate((simu_ensemble_out,simu_ensemble[0:obs_temp.ntime-1,:]),axis=0)
            obs_out = np.concatenate((obs_out,np.concatenate((obs_temp.time[0:obs_temp.ntime-1].reshape(-1,1),obs_temp.value[0:obs_temp.ntime-1,:]),axis=1)),axis=0)
            
        state_vector_out = np.concatenate((state_vector_out,state_vector),axis=0)
        if with_head:
            np.savetxt('./pflotran_results/state_vector_out_with_head.txt',state_vector_out)
            np.savetxt('./pflotran_results/simu_ensemble_out_with_head.txt',simu_ensemble_out)
            np.savetxt('./pflotran_results/obs_out_with_head.txt',obs_out)
        else:
            np.savetxt('./pflotran_results/state_vector_out_without_head.txt',state_vector_out) 
            np.savetxt('./pflotran_results/simu_ensemble_out_without_head.txt',simu_ensemble_out)
            np.savetxt('./pflotran_results/obs_out_without_head.txt',obs_out)
            
        mod.state_vector = inflate_state_vector(state_vector,infl)
        ibatch = ibatch+1
    return state_vector_out

def inflate_state_vector(state_vector,infl):
    state_mean = np.mean(state_vector,axis=1,keepdims=True)
    state_diff = state_vector-state_mean
    state_vector = state_mean+state_diff*infl

    return state_vector
