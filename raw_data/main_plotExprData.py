
import matplotlib.pyplot as plt
import numpy as np
import metadata_monks as mm
import my_plots as mp

import warnings
warnings.filterwarnings("ignore")

monkey_ID = 2

import sys, os
folder_input_data = os.getcwd() + '/abunds/'

#[Mss_QVI, sample_ratio, months] = rs.get_monkey_parameters(monkey_ID=monkey_ID,SKIP_MONTH2=False)
grans, grans_EGFP, S_QVIs, months = mm.get_meta_data(monkey_ID)
print 'S_QVIs:', S_QVIs
#data_expr = rs.get_expr_data(monkey_ID=monkey_ID,SKIP_MONTH2=False)
data_expr = np.loadtxt(folder_input_data+'data_expr_'+str(monkey_ID)+'.txt')
row, col = data_expr.shape

## how many clones show up in each sample? ##
#for j in range(col):
#    print sum(data_expr[:,j] > 0)
#############################################

TO_PLOT_PAIRWISE = False
if TO_PLOT_PAIRWISE:
    data_expr = data_expr[:,1:]
    row, col = data_expr.shape
    all_pairs = []
    pair_stats = []
    for i in range(row):
        for j in range(col-1):
            if data_expr[i,j] == 0:
                continue
            all_pairs.append([data_expr[i,j],data_expr[i,j+1]])
    all_pairs = np.array(all_pairs)
    plt.scatter(np.log(all_pairs[:,0]), np.log(all_pairs[:,1]))
    plt.show()
    quit()

SORT_FREQS = True
if SORT_FREQS:
    ''
    avgs_all = []
    for i in range(row):
        curr_avg = data_expr[i,:].mean()
        avgs_all.append(curr_avg)
    inds_sorted = np.argsort(avgs_all)
    inds_sorted[:] = inds_sorted[::-1] # small to big
    data_expr = data_expr[inds_sorted,:]

plot_diversity = True
if plot_diversity:
    cell_nums = np.zeros(data_expr.shape)
    richness = []
    simpsons = []
    shannons = []
    for j in range(data_expr.shape[1]):
        cur_simp = sum(data_expr[:, j] ** 2)
        simpsons.append(cur_simp*50.)

        ha = data_expr[:, j] * np.log(data_expr[:, j])
        cell_nums[:, j] = data_expr[:, j] * S_QVIs[j]
        for i in range(data_expr.shape[0]):
            cell_nums[i, j] = int(round(cell_nums[i, j]))
            if ha[i] != ha[i]:
                ha[i] = 0
        cur_shan = -sum(ha)
        shannons.append(cur_shan)

        cur_rich = sum(cell_nums[:, j] > 0)
        richness.append(cur_rich/50.)
    fig = plt.figure(figsize=(8, 6))

    plt.plot(months, shannons, 'k', linewidth=1.5, label='Shannon')
    plt.plot(months, richness, 'b', linewidth=1.5, label='0.02 * Richness')
    plt.plot(months, simpsons, 'g', linewidth=1.5, label='50 * Simpson')

    mp.config_plot('$t_j$ (in month)', 'diversity index', legend_loc=(0.9, 0.4))
    mp.use_sci_nota(False, False)
    plt.xlim(months[0], months[-1])
    plt.ylim(0, 5)
    plt.show()
    quit()

cumu_fhi = []
for j in range(col):
    curr_cumu = np.cumsum(data_expr[:,j])
    cumu_fhi.append(curr_cumu)
cumu_fhi = np.array(cumu_fhi)

####### Plot Abundances #######
TO_PLOT_ABUND = True
if TO_PLOT_ABUND:
    avg_par = []
    for i in range(row):
        avg_par.append(data_expr[i,1:].mean())
    avg_par = np.array(avg_par)
    large_marks_pre = avg_par > 0.012 # thicker curve
    small_marks_pre = avg_par < 0.0005 # thinner curve
    #inter_marks_pre = (avg_par >= small_marks_pre) and (avg_par <= large_marks_pre)

    fig = plt.figure(figsize=(8,6))
    num_abund_clones = cumu_fhi.shape[1]
    linewidths = np.array([.5]*row) # pre-set to intermediate clones setting
    linewidths[large_marks_pre] = np.array([1.]*sum(large_marks_pre))
    linewidths[small_marks_pre] = np.array([.1]*sum(small_marks_pre))

    for j in range(num_abund_clones):
        plt.plot(months, cumu_fhi[:,j], linewidth=linewidths[j])
    mp.config_plot('$t_j$ (in month)', '$\hat{f}_{i}(t_j)$ (stacked)')
    mp.use_sci_nota(False, False)
    plt.xlim(months[0], months[-1])
    plt.ylim(0,1)

    plt.show()
    quit()

grans, grans_EGFP, S_QVIs, months = mm.get_meta_data(monkey_ID)

TO_PLOT_TOTAL = False
if TO_PLOT_TOTAL:
    rnd_size = 180

    fig, ax = plt.subplots(figsize=(8,6))
    #plt.plot(months, tot_blood, color ='b', linestyle='-',
    #         label='Total Grans+PBMC', linewidth=1)
    #plt.scatter(months, tot_blood, edgecolor='b', facecolors='b', marker='^', s=rnd_size)

    plt.plot(months, grans, color='b', linestyle='-',
             label='Total Grans', linewidth=1.5)
    plt.scatter(months, grans, edgecolor='b', facecolors='b', marker='^', s=rnd_size)

    plt.plot(months, grans_EGFP, color='g', linestyle='-',
             label='EGFP+ Grans', linewidth=2.5)
    plt.scatter(months, grans_EGFP, edgecolor='g', facecolors='g', marker='s', s=rnd_size)

    plt.plot(months, S_QVIs, color='k', linestyle='-',
             label='EGFP+ Grans', linewidth=2.5)
    plt.scatter(months, S_QVIs, edgecolor='k', facecolors='k', marker='o', s=rnd_size)

    ax.set_yscale('log')
    mp.config_plot('$t_j$ (in month)', 'sample sizes') #, legend_loc=(.9,.5))
    ax.set_xlim(months[0], months[-1])
    ax.set_ylim(10**3, 4*10**6)
    plt.show()

TO_PLOT_SAMPLING_VAR = True
if TO_PLOT_SAMPLING_VAR:
    row, col = data_expr.shape
    avgs = np.mean(data_expr, axis=1)
    stds = np.std(data_expr, axis=1)

    plt.scatter(np.log(avgs),np.log(stds),edgecolors='b',facecolors='none')


    S = S_QVIs.mean() #7.7*10.**3
    #print S
    bnd_large = 0.1
    bnd_small = 0.0001
    p_s = np.array([np.exp(x) for x in np.linspace(np.log(bnd_small),np.log(bnd_large),num=100)])

    data_theo = []
    for p in p_s:
        one_theo_clone = np.random.binomial(S, p, size=8)
        data_theo.append(one_theo_clone/S)
    data_theo = np.array(data_theo)

    row, col = data_theo.shape

    avgs = []
    stds = []
    j0, j1 = 5, col
    for i in range(row):
        curr_mean = data_theo[i,:].mean()
        if curr_mean < bnd_small:
            continue

        avgs.append(data_theo[i,j0:j1].mean())
        stds.append(data_theo[i,j0:j1].std(ddof=1))
    avgs = np.array(avgs)
    stds = np.array(stds)

    ax = plt.subplot(111)
    plt.scatter(np.log(avgs), (stds), label='simulated \n sampling',linewidth=1.5,
                edgecolors='g',facecolors='none', marker='^', s=120)
    plt.xlim([np.log(bnd_small),np.log(bnd_large)])
    plt.show()