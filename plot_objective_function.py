import numpy as np
import matplotlib.pyplot as plt

def main():
    path_to_files = '../RunEkoProdDistLoc'
    num_iter = 4

    mm_dist = []
    for iter in range(num_iter):
        mm_dist.append(np.load(str(path_to_files) + '/debug_analysis_step_{}.npz'.format(iter))['data_misfit'])

    # correlation based
    path_to_files = '../RunEkoProdDistLoc_CrossValidation'

    mm_corr = []
    for iter in range(num_iter):
        mm_corr.append(np.load(str(path_to_files) + '/debug_analysis_step_{}.npz'.format(iter))['data_misfit'])

    plt.figure(); plt.plot(mm_corr,'bo-', label='Corr_loc');plt.plot(mm_dist,'ko-', label='Dist_loc');
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 16)
    plt.xlabel('Iteration no.', size=20)
    plt.ylabel('Data mismatch', size=20)
    # plt.title('Objective function', size=28)
    plt.show()





main()