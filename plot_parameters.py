import numpy as np
import matplotlib.pyplot as plt
import pickle

def plot_para():

    path_to_files = '../RunEkoProdDistLoc'
    save_figure = False

    f_dim= (22,155,128)
    actnum = np.load(str(path_to_files)+'/actnum.npz')['actnum']

    layer_no = 10

    # post = np.load(str(path_to_files)+'/final_ens4_adaloc.npz')['permx']

    # ##### Load True permx
    # true_perm = []
    # with open(str(path_to_files)+'/True/permx.dat','r') as f:
    #     lines = f.readlines()
    #     for el in lines:
    #         true_perm.append(el)
    #
    #
    # permx_true = []
    # for el in true_perm[1:-1]:
    #     for val in el.split():
    #          permx_true.append(np.float(val))
    #
    # permx_true = np.log(np.array(permx_true))
    # permx_true[~actnum]= np.nan
    # permx_true = np.reshape(permx_true, f_dim)

    ##### Load Prior permx and std

    prior = np.load(str(path_to_files) + '/prior.npz')['permx']
    permx_prior = np.zeros(f_dim)
    permx_prior_std = np.zeros(f_dim)
    permx_prior[:] = np.nan
    permx_prior_std[:] = np.nan
    permx_prior[actnum.reshape(f_dim)] = prior.mean(1)
    permx_prior_std[actnum.reshape(f_dim)] = prior.std(axis=1, ddof=1)

    ##### Load Post permx and std

    post = np.load(str(path_to_files) + '/final_distLoc.npz')['permx']
    permx_post = np.zeros(f_dim)
    permx_post_std = np.zeros(f_dim)
    permx_post[:] = np.nan
    permx_post_std[:] = np.nan
    permx_post[actnum.reshape(f_dim)] = post.mean(1)
    permx_post_std[actnum.reshape(f_dim)] = post.std(axis=1, ddof=1)

    # ####### Plotting
    Colorbar_max = 7.0
    Colorbar_min = 0.0

    plt.figure();plt.imshow(permx_prior[layer_no,:,:], cmap='jet', vmax= Colorbar_max, vmin=Colorbar_min);plt.colorbar(); plt.title('Prior permx (log) at layer ' + str(layer_no + 1), size=20 )
    plt.xticks(fontsize = 14); plt.yticks(fontsize = 14)
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/Prior_Permx(log)_Layer_' + str(layer_no + 1))

    plt.figure();plt.imshow(permx_post[layer_no,:,:], cmap='jet', vmax= Colorbar_max, vmin=Colorbar_min);plt.colorbar(); plt.title('Post permx (log) at layer ' + str(layer_no + 1) + ', Dist_loc', size=20 )
    plt.xticks(fontsize = 14); plt.yticks(fontsize = 14)
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/Post_Permx(log)_Layer_' + str(layer_no + 1))

    # plt.figure();plt.imshow(permx_true[layer_no,:,:], cmap='jet', vmax= Colorbar_max, vmin=Colorbar_min);plt.colorbar(); plt.title('True permx (log) at layer ' + str(layer_no + 1), size=28 )
    # plt.xticks(fontsize = 16); plt.yticks(fontsize = 16)

    # plt.show()


     # Plot stds
    plt.figure(); plt.imshow(permx_prior_std[layer_no, :, :], cmap='jet')
    plt.xticks(fontsize = 14); plt.yticks(fontsize = 14)
    plt.colorbar(); plt.title('Prior std permx (log) at layer ' + str(layer_no + 1), size=20 )
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/Prior_Std_Permx(log)_Layer_' + str(layer_no + 1))

    plt.figure(); plt.imshow(permx_post_std[layer_no, :, :], cmap='jet')
    plt.xticks(fontsize = 14); plt.yticks(fontsize = 14)
    plt.colorbar(); plt.title('Post std permx (log) at layer ' + str(layer_no + 1) + ', Dist_loc', size=20 )
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/Post_Std_Permx(log)_Layer_' + str(layer_no + 1))

    plt.show()


plot_para()