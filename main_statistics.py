

import numpy as np
import os
import utils_simulations as rs
import re

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

from main_simulations import L_s, r2_s

simulated_folderpath = 'simulated_clones/'

res_2D_folderpath = 'figs/Fig_2D_surf/'
res_1D_folderpath = 'figs/Fig_1D_line/'
res_pt_folderpath = 'figs/Fig_pt_MSE/'
res_pt_avgstd_folderpath = 'figs/Fig_pt_MSE_avgstd/'

def plot_2D_surface(use_saved_res = True):

    if not os.path.exists(res_2D_folderpath):
        os.mkdir(res_2D_folderpath)

    if not use_saved_res:
        X, Y = np.meshgrid(r2_s, L_s)
        Z = np.zeros((len(r2_s),len(L_s)))

        [Mss_EGFP, sample_ratio, months] = rs.get_monkey_parameters(monkey_ID=3)

        data_expr = rs.get_expr_data(monkey_ID=3)
        fz_expr_avg, fz_expr_std = rs.get_fz_stats(data_expr)

        repeat_num = 200

        all_dists = []
        for i in range(len(r2_s)):
            r2 = r2_s[i]
            for j in range(len(L_s)):
                L = L_s[j]

                filename = simulated_folderpath + '[r2,L]='+"[%.1f,%.1f]"%(r2,L)
                print filename

                mi_data = np.loadtxt(filename)
                mi_row, mi_col = mi_data.shape
                sample_days_inds=range(3*30,mi_col,60)

                sample_big = rs.get_samples(mi_s=mi_data,
                                              sample_days_inds=sample_days_inds,
                                              sample_ratio=sample_ratio)

                curr_dists = []
                for k in range(repeat_num):
                    samples_theo = rs.resample(sample_big, 8)
                    fz_theo_avg, fz_theo_std = rs.get_fz_stats(samples_theo)
                    curr_dists.append(np.square(rs.my_euclid(fz_theo_avg[1:],fz_expr_avg[1:])))

                all_dists.append(curr_dists)
                curr_dists = np.array(curr_dists)
                Z[i,j] = curr_dists.mean()
        all_dists = np.array(all_dists)

        np.savetxt(res_2D_folderpath+'X', X)
        np.savetxt(res_2D_folderpath+'Y', Y)
        np.savetxt(res_2D_folderpath+'Z', Z)

        np.savetxt(res_2D_folderpath+'all_dists', all_dists)
    else:
        X = np.loadtxt(res_2D_folderpath+'X')
        Y = np.loadtxt(res_2D_folderpath+'Y')
        Z = np.loadtxt(res_2D_folderpath+'Z')

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib as mpl
    mpl.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib import cm


    #thres = 0.005
    #Z[Z>thres] = thres


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X,Y,np.transpose(Z),
                        rstride=1, cstride=1, alpha=0.5, linewidth=0, cmap=cm.jet,antialiased=True)
    fig.colorbar(surf,shrink=.5,aspect=5)
    fig.subplots_adjust(top=1, bottom=0, left=0, right=1)

    ax.set_xlabel('$r_n$', fontsize=20)
    ax.set_ylabel('$L$', fontsize=20)
    ax.set_zlabel('energy', fontsize=25)

    ax.w_zaxis.line.set_lw(0.)
    ax.set_zticks([])

    #ax.set_axis_off()
    ax.view_init(azim=0, elev=90)
    plt.savefig(res_2D_folderpath + '2D_surf.png')



def plot_1D_line(use_saved_files = True):

    r2 = 2.5

    from main_simulations import L_s, res_folderpath

    [Mss_EGFP, sample_ratio, months] = rs.get_monkey_parameters(monkey_ID=3)

    change_Mss = False
    if change_Mss:
        Mss_EGFP /= 5.
        sample_ratio *= 5

    data_expr = rs.get_expr_data(monkey_ID=3)
    fz_expr_avg, fz_expr_std = rs.get_fz_stats(data_expr)



    if not use_saved_files:
        repeat_num = 200

        all_dists = []
        all_avgs = []
        all_stds = []

        for j in range(len(L_s)):
            L = L_s[j]

            filename = '[r2,L]=' + "[%.1f,%.1f]" % (r2, L)
            print filename

            mi_data = np.loadtxt(res_folderpath + filename)
            mi_row, mi_col = mi_data.shape
            sample_days_inds = range(3 * 30, mi_col, 60)

            sample_big = rs.get_samples(mi_s=mi_data,
                                        sample_days_inds=sample_days_inds,
                                        sample_ratio=sample_ratio)

            curr_dists = []
            for k in range(repeat_num):
                samples_theo = rs.resample(sample_big, 8)
                fz_theo_avg, fz_theo_std = rs.get_fz_stats(samples_theo)
                curr_dists.append(np.square(rs.my_euclid(fz_theo_avg[1:], fz_expr_avg[1:])))

            all_dists.append(curr_dists)
            curr_dists = np.array(curr_dists)

            all_avgs.append(curr_dists.mean())
            all_stds.append(curr_dists.std())

        if not change_Mss:
            np.savetxt(res_1D_folderpath+'all_avgs_fix_r2', all_avgs)
            np.savetxt(res_1D_folderpath+'all_stds_fix_r2', all_stds)
        else:
            np.savetxt(res_1D_folderpath+'all_avgs_fix_r2', all_avgs)
            np.savetxt(res_1D_folderpath+'all_stds_fix_r2', all_stds)
    else:
        if not change_Mss:
            all_avgs = np.loadtxt(res_1D_folderpath+'all_avgs_fix_r2')
            all_stds = np.loadtxt(res_1D_folderpath+'all_stds_fix_r2')
        else:
            all_avgs = np.loadtxt(res_1D_folderpath+'all_avgs_fix_r2')
            all_stds = np.loadtxt(res_1D_folderpath+'all_stds_fix_r2')

    import my_funcs as mf
    plt.errorbar(L_s[:-4], all_avgs[:-4], yerr=all_stds[:-4],
                 color='k', linewidth=2)
    mf.config_plot(r'$L_{\rm e}$', 'MSE')
    mf.use_sci_nota('x', usesci=False)
    mf.use_sci_nota('y', usesci=True)

    print L_s[:-4][np.argmin(all_avgs[:-4])]

    plt.ylim([0, 0.00032])  # plt.ylim([0,0.00025]) #
    plt.xlim([18.8, 26.2])
    # plt.yticks(np.linspace(0.005,0.025,num=5).tolist())
    # plt.tight_layout()
    plt.savefig(res_1D_folderpath+'fig7a.png')  # , format='svg', dpi=600
    # plt.show()

    ind_min = np.argmin(all_avgs)
    avg_min = all_avgs[ind_min]
    std_min = all_stds[ind_min]
    L_min = L_s[ind_min]

    print avg_min, std_min, L_min


def plot_pt_MSE(use_saved_files = True):
    L = 24.0  # 23.4
    r2 = 2.5

    [Mss_EGFP, sample_ratio, months] = rs.get_monkey_parameters(monkey_ID=3)

    data_expr = rs.get_expr_data(monkey_ID=3)
    row, col = data_expr.shape

    fz_expr_avg, fz_expr_std = rs.get_fz_stats(data_expr)


    if not os.path.exists(res_pt_folderpath):
        os.mkdir(res_pt_folderpath)

    if not use_saved_files:

        all_dists = []
        all_avgs = []
        all_stds = []

        filename = 'simulated_clones/' + '[r2,L]=' + "[%.1f,%.1f]" % (r2, L)
        print filename

        mi_data = np.loadtxt(filename)
        mi_row, mi_col = mi_data.shape
        sample_days_inds = range(3 * 30, mi_col, 60)

        sample_big = rs.get_samples(mi_s=mi_data,
                                    sample_days_inds=sample_days_inds,
                                    sample_ratio=sample_ratio)

        repeat_num = 200

        all_stats = []
        for k in range(repeat_num):
            samples_theo = rs.resample(sample_big, col)
            all_stats += rs.two_sample_test(samples_theo, data_expr)[:-1]
        all_stats = np.array(all_stats)
        print sum(all_stats >= 0.05) / float(len(all_stats))

        fz_theo_avg, fz_theo_std = rs.get_fz_stats(samples_theo)
        np.savetxt(res_pt_folderpath+'fz_theo_avg.txt', fz_theo_avg)
        np.savetxt(res_pt_folderpath+'fz_theo_std.txt', fz_theo_std)
    else:
        fz_theo_avg = np.loadtxt(res_pt_folderpath+'fz_theo_avg.txt')
        fz_theo_std = np.loadtxt(res_pt_folderpath+'fz_theo_std.txt')

    import my_funcs as mf
    fig = plt.figure(figsize=(8, 6))
    plt.errorbar(range(1, 8), fz_expr_avg[1:], yerr=fz_expr_std[1:],
                 color='k', linewidth=2)
    eb_theo = plt.errorbar(range(1, 8), fz_theo_avg[1:], yerr=fz_theo_std[1:],
                           color='k', linewidth=2, linestyle='--')
    eb_theo[-1][0].set_linestyle('--')
    mf.config_plot(r'$z$', r'$\qquad Y_{z}(L_{\rm e}), ~~~~~~\hat{Y}_{z}$')
    mf.use_sci_nota('x', usesci=False)
    mf.use_sci_nota('y', usesci=True)
    plt.xlim([0.5, 7.5])
    plt.ylim([-0.0005, 0.007])
    # plt.yticks(np.linspace(0.005,0.025,num=5).tolist())
    # plt.tight_layout()

    plt.savefig(res_pt_folderpath+'fig7b.png') #, format='svg', dpi=600


def plot_pt_MSE_avgstd(use_saved_files = False):
    L = 24.0 #23.9
    r2 = 2.5

    [Mss_EGFP, sample_ratio, months] = rs.get_monkey_parameters(monkey_ID=3)

    data_expr = rs.get_expr_data(monkey_ID=3)
    row, col = data_expr.shape

    def get_avgstds(data):
        return np.mean(data, axis=1), np.std(data, axis=1)

    fz_expr_avg, fz_expr_std = get_avgstds(data_expr)

    if not os.path.exists(res_pt_avgstd_folderpath):
        os.mkdir(res_pt_avgstd_folderpath)


    if use_saved_files and os.path.isfile(res_pt_avgstd_folderpath+'samples_theo_' + str(L) + '.txt'):
        samples_theo = np.loadtxt(res_pt_avgstd_folderpath+'samples_theo_' + str(L) + '.txt')
    else:

        filename = simulated_folderpath + '[r2,L]=' + "[%.1f,%.1f]" % (r2, L)
        print filename

        mi_data = np.loadtxt(filename)
        mi_row, mi_col = mi_data.shape
        sample_days_inds = range(3 * 30, mi_col, 60)

        sample_big = rs.get_samples(mi_s=mi_data,
                                    sample_days_inds=sample_days_inds,
                                    sample_ratio=sample_ratio)

        repeat_num = 1

        all_stats = []
        for k in range(repeat_num):
            samples_theo = rs.resample(sample_big, col)
            all_stats += rs.two_sample_test(samples_theo, data_expr)[:-1]
        all_stats = np.array(all_stats)
        print sum(all_stats >= 0.05) / float(len(all_stats))

        np.savetxt(res_pt_avgstd_folderpath+'samples_theo_' + str(L) + '.txt', samples_theo)

    TO_PLOT_ABUND = False
    if TO_PLOT_ABUND:
        import my_funcs as mf

        cumu_fhi = []
        for j in range(col):
            curr_cumu = np.cumsum(samples_theo[:, j])
            cumu_fhi.append(curr_cumu)
        cumu_fhi = np.array(cumu_fhi)
        PRE_CLASSIFY_CLONES = True
        avg_par = []
        for i in range(row):
            avg_par.append(data_expr[i, 1:].mean())
        avg_par = np.array(avg_par)
        large_marks_pre = avg_par > 0.012
        small_marks_pre = avg_par < 0.0005
        # inter_marks_pre = (avg_par >= small_marks_pre) and (avg_par <= large_marks_pre)

        fig = plt.figure(figsize=(8, 6))
        num_abund_clones = cumu_fhi.shape[1]
        linewidths = np.array([.5] * row)  # pre-set to intermediate clones setting
        linewidths[large_marks_pre] = np.array([1.] * sum(large_marks_pre))
        linewidths[small_marks_pre] = np.array([.1] * sum(small_marks_pre))

        # [exp(x) for x in avgs_all] #[exp(x) for x in np.linspace(log(1.),log(.1),num_abund_clones)]
        # print linewidths
        for j in range(num_abund_clones):
            plt.plot(mf.months, cumu_fhi[:, j], linewidth=linewidths[j])
        # plt.plot(mf.months, cumu_fhi, linewidth=.7)
        mf.config_plot('$t_j$ (in month)', '$\hat{f}_{i}(t_j)$ (stacked)')
        mf.use_sci_nota('x', usesci=False)
        mf.use_sci_nota('y', usesci=False)
        plt.xlim(mf.months[0], mf.months[-1])
        plt.ylim(0, 1)

        # plt.axis('off')
        plt.show()

    def is_large_clone(data):
        res = []
        row, col = data.shape
        for i in range(row):
            res.append(sum(data[i, :] == 0) == 0)  # if clone i has some absence
        return np.array(res)

    fz_theo_avg, fz_theo_std = get_avgstds(samples_theo)
    large_mrks = is_large_clone(samples_theo)
    other_mrks = [not x for x in large_mrks]

    import my_funcs as mf
    plt.scatter(np.log(fz_expr_avg), fz_expr_std, c='b', label='experiment')

    ax = plt.subplot(111)
    for_ppt = False

    if for_ppt:
        curr_label = 'model + sampling'
    else:
        curr_label = 'model + \n sampling'  # 'simulated \n sampling'
    plt.scatter(np.log(fz_theo_avg), fz_theo_std, label=curr_label, linewidth=1.5,
                edgecolors='g', facecolors='none', marker='^', s=120)
    # plt.scatter(np.log(fz_theo_avg[large_mrks]), fz_theo_std[large_mrks], label='simulated \n sampling',linewidth=1.5,
    #                 edgecolors='r',facecolors='none', marker='^', s=120)
    plt.xlim([-8.5, -2.1])
    plt.ylim([-0.005, 0.05])  # plt.ylim([-0.001,0.01]) #

    # mf.config_plot(r'$\ln y_i$',r'$\sigma_i$',legend_loc=(.65,.8))

    if not for_ppt:
        mf.config_plot('', '', legend_loc=(.65, .8))
    else:
        mf.config_plot(r'ln($y_i$)', r'$\sigma_i$', legend_loc=(.7, .8))
    mf.use_sci_nota('x', False)
    mf.use_sci_nota('y', True)
    # plt.show()
    plt.savefig(res_pt_avgstd_folderpath+'fig9b.png') #, format='svg', dpi=600


if __name__ == '__main__':
    # plot_pt_MSE(use_saved_files = True)
    # plot_1D_line(use_saved_files = True)
    # plot_2D_surface(use_saved_res = True)
    plot_pt_MSE_avgstd(use_saved_files = False)