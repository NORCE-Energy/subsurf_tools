import numpy as np
import matplotlib.pyplot as plt
import pickle


def new_plot():
    path_to_files = '../RunEkofisk_Production'

    n = 200
    ne=4



    obs = np.load(str(path_to_files)+'/obs_var.npz')['obs']
    var = np.load(str(path_to_files) + '/obs_var.npz')['var']

    tot_key = [el for el in obs[0].keys()]


    with open(str(path_to_files)+ '/sim_results.p','rb') as f:
        pred_data = pickle.load(f)


    ##### Oil production rate

    my_data = tot_key[n]
    t1, t2 = my_data.split()

    data1 = [pred_data[i][my_data] if pred_data[i][my_data] is not None else np.ones((1, ne)) * np.nan for i in range(919)]
    data1 = np.concatenate(data1,axis=0)
    obs_data = [obs[i][my_data] if obs[i][my_data] is not None else np.nan for i in range(919)]
    stddev = [np.sqrt(var[i][my_data]) if obs[i][my_data] is not None else np.nan for i in range(919)]

    plt.figure();
    plt.plot(data1, 'k');
    plt.plot(obs_data, 'r');
    plt.plot([obs_data[el] + stddev[el] for el in range(919)], ':r');
    plt.plot([obs_data[el] - stddev[el] for el in range(919)], ':r');
    plt.title(str(t1) + ' forcast, at Well: ' + str(t2))


    ##### Water production rate

    my_data = tot_key[n + 312]
    t1, t2 = my_data.split()

    data1 = [pred_data[i][my_data] if pred_data[i][my_data] is not None else np.ones((1, ne)) * np.nan for i in range(919)]
    data1 = np.concatenate(data1,axis=0)
    obs_data = [obs[i][my_data] if obs[i][my_data] is not None else np.nan for i in range(919)]
    stddev = [np.sqrt(var[i][my_data]) if obs[i][my_data] is not None else np.nan for i in range(919)]

    plt.figure();
    plt.plot(data1, 'k');
    plt.plot(obs_data, 'r');
    plt.plot([obs_data[el] + stddev[el] for el in range(919)], ':r');
    plt.plot([obs_data[el] - stddev[el] for el in range(919)], ':r');
    plt.title(str(t1) + ' forcast, at Well: ' + str(t2))

    ##### Gas production rate

    my_data = tot_key[n + 624]
    t1, t2 = my_data.split()

    data1 = [pred_data[i][my_data] if pred_data[i][my_data] is not None else np.ones((1, ne)) * np.nan for i in range(919)]
    data1 = np.concatenate(data1,axis=0)
    obs_data = [obs[i][my_data] if obs[i][my_data] is not None else np.nan for i in range(919)]
    stddev = [np.sqrt(var[i][my_data]) if obs[i][my_data] is not None else np.nan for i in range(919)]

    plt.figure();
    plt.plot(data1, 'k');
    plt.plot(obs_data, 'r');
    plt.plot([obs_data[el] + stddev[el] for el in range(919)], ':r');
    plt.plot([obs_data[el] - stddev[el] for el in range(919)], ':r');
    plt.title(str(t1) + ' forcast, at Well: ' + str(t2))



    plt.show()

def main():

    path_to_files = '../RunEkoProdDistLoc'  ####  RunEkoProdCorrLoc
    save_figure = True  ## Use True  for saving the figures

    n = 300 ## Select a well
    ne = 100
    ###########

    obs = np.load(str(path_to_files)+'/obs_var.npz')['obs']
    test1 = np.load(str(path_to_files)+'/debug_analysis_step_0.npz')
    test2 = np.load(str(path_to_files)+'/debug_analysis_step_4.npz')

    tot_key = [el for el in obs[0].keys()]
    my_data = tot_key[n]
    t1,t2=my_data.split()

    ##### Oil production rate

    n_d_obs_o = [obs[i][my_data] for i in range(919)]


    pred1 = test1['pred_data']
    data1 = [pred1[i][my_data] if pred1[i][my_data] is not None else np.ones((1,ne))*np.nan for i in range(919)]
    n_d1 = np.concatenate(data1,axis=0)

    pred2 = test2['pred_data']
    data2 = [pred2[i][my_data] if pred2[i][my_data] is not None else np.ones((1,ne))*np.nan for i in range(919)]
    n_d2 = np.concatenate(data2,axis=0)

    plt.figure();plt.plot(n_d1,'k');plt.plot(n_d_obs_o,'r'); plt.title(str(t1) + ' forcast for iteration = 0, at Well: ' + str(t2))
    ylim = plt.gca().get_ylim()
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/' + str(t2[0:2]) + str(t2[3:]) + '_' + str(t1) + '_0')
    plt.figure();plt.plot(n_d2,'k');plt.plot(n_d_obs_o,'r'); plt.title(str(t1) + ' forcast for iteration = 4, at Well: ' + str(t2))
    plt.gca().set_ylim(ylim)
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/' + str(t2[0:2]) + str(t2[3:]) + '_' + str(t1) + '_4')

     ##### Water production rate
    my_data = tot_key[n + 312]
    t1,t2=my_data.split()

    n_d_obs_w = [obs[i][my_data] for i in range(919)]

    pred1 = test1['pred_data']
    data1_w = [pred1[i][my_data] if pred1[i][my_data] is not None else np.ones((1,ne))*np.nan for i in range(919)]
    n_d1 = np.concatenate(data1_w,axis=0)

    pred2 = test2['pred_data']
    data2 = [pred2[i][my_data] if pred2[i][my_data] is not None else np.ones((1,ne))*np.nan for i in range(919)]
    n_d2 = np.concatenate(data2,axis=0)

    plt.figure();plt.plot(n_d1,'k');plt.plot(n_d_obs_w,'r'); plt.title(str(t1) + ' forcast for iteration = 0, at Well: ' + str(t2))
    ylim = plt.gca().get_ylim()
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/' + str(t2[0:2]) + str(t2[3:]) + '_' + str(t1) + '_0')
    plt.figure();plt.plot(n_d2,'k');plt.plot(n_d_obs_w,'r'); plt.title(str(t1) + ' forcast for iteration = 4, at Well: ' + str(t2))
    plt.gca().set_ylim(ylim)
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/' + str(t2[0:2]) + str(t2[3:]) + '_' + str(t1) + '_4')
    # plt.show()

     ##### Gas production rate
    my_data = tot_key[n + 624]
    t1,t2=my_data.split()

    n_d_obs_g = [obs[i][my_data] for i in range(919)]

    pred1 = test1['pred_data']
    data1 = [pred1[i][my_data] if pred1[i][my_data] is not None else np.ones((1,ne))*np.nan for i in range(919)]
    n_d1 = np.concatenate(data1,axis=0)

    pred2 = test2['pred_data']
    data2 = [pred2[i][my_data] if pred2[i][my_data] is not None else np.ones((1,ne))*np.nan for i in range(919)]
    n_d2 = np.concatenate(data2,axis=0)

    plt.figure();plt.plot(n_d1,'k');plt.plot(n_d_obs_g,'r');plt.title(str(t1) + ' forcast for iteration = 0, at Well: ' + str(t2))
    ylim = plt.gca().get_ylim()
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/' + str(t2[0:2]) + str(t2[3:]) + '_' + str(t1) + '_0')
    #plt.figure();plt.plot(n_d_obs_o,'k');plt.plot(n_d_obs_g,':m');plt.plot(n_d_obs_w,'--r');
    plt.figure();plt.plot(n_d2,'k');plt.plot(n_d_obs_g,'r'); plt.title(str(t1) + ' forcast for iteration = 4, at Well: ' + str(t2))
    plt.gca().set_ylim(ylim)
    if save_figure is True: plt.savefig(str(path_to_files) +'/Figures/' + str(t2[0:2]) + str(t2[3:]) + '_' + str(t1) + '_4')
    #plt.show()

    ############
    plt.show()

main()
# new_plot()