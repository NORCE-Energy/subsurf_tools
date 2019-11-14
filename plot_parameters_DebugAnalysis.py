import numpy as np
import matplotlib.pyplot as plt
import pickle

def plot_para():

    path_to_files = '../RunEkoProdDistLoc_CrossValidation'
    save_figure = True  ## Use True  for saving the figures

    f_dim= (22,155,128)
    actnum = np.load(str(path_to_files)+'/actnum.npz')['actnum']

    layer_no = 9

    # Load debug steps
    num_iter = 4
    for iter in range(num_iter):
        tmp_post = np.zeros(f_dim)
        tmp_post[:] = np.nan
        tmp_post[actnum.reshape(f_dim)] = \
            np.load(str(path_to_files) + '/debug_analysis_step_{}.npz'.format(iter))['state'][()]['permx'].std(axis=1, ddof=1)

        plt.figure();
        plt.imshow(tmp_post[layer_no, :, :], cmap='jet');
        plt.colorbar();
        plt.title('Iter {} std. permx (log) at layer '.format(iter) + str(layer_no))

    ##### Load Prior permx

    prior= np.load(str(path_to_files)+'/prior.npz')['permx']
    permx_prior = np.zeros(f_dim)
    permx_prior_std = np.zeros(f_dim)
    permx_prior[:]=np.nan
    permx_prior_std[:]=np.nan
    permx_prior[actnum.reshape(f_dim)] = prior.mean(1)
    permx_prior_std[actnum.reshape(f_dim)] = prior.std(axis=1,ddof=1)

    ##### Load Post permx

    post = np.load(str(path_to_files) + '/final_adaLoc.npz')['permx']
    permx_post = np.zeros(f_dim)
    permx_post_std = np.zeros(f_dim)
    permx_post[:] = np.nan
    permx_post_std[:] = np.nan
    permx_post[actnum.reshape(f_dim)] = post.mean(1)
    permx_post_std[actnum.reshape(f_dim)] = post.std(axis=1, ddof=1)


     ####### Plotting
    Colorbar_max = 8.0
    Colorbar_min = 0.0

    plt.figure();plt.imshow(permx_prior[layer_no,:,:], cmap='jet', vmax= Colorbar_max, vmin=Colorbar_min);plt.colorbar(); plt.title('Prior permx (log) at layer ' + str(layer_no))
    plt.figure();plt.imshow(permx_post[layer_no,:,:], cmap='jet', vmax= Colorbar_max, vmin=Colorbar_min);plt.colorbar(); plt.title('Post permx (log) at layer ' + str(layer_no))

    # plt.figure();plt.imshow(permx_post[layer_no,:,:] - permx_prior[layer_no,:,:] , cmap='jet');plt.colorbar(); plt.title('Diff permx (log) at layer ' + str(layer_no))
    # plt.figure();plt.imshow(np.exp(permx_post[layer_no,:,:]) - np.exp(permx_prior[layer_no,:,:]) , cmap='jet');plt.colorbar(); plt.title('Diff permx at layer ' + str(layer_no))

    # std
    plt.figure();
    plt.imshow(permx_prior_std[layer_no, :, :], cmap='jet');
    plt.colorbar();
    plt.title('Prior std permx (log) at layer ' + str(layer_no))

    plt.figure();
    plt.imshow(permx_post_std[layer_no, :, :], cmap='jet');
    plt.colorbar();
    plt.title('Post std permx (log) at layer ' + str(layer_no))

    plt.show()


plot_para()