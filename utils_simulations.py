
import numpy as np
import random as rd
from scipy.special import gammainc
from numpy.random import binomial
#import pandas as pd
from scipy.spatial.distance import euclidean
import math

data_folderpath = 'raw_data/'

def get_expr_data(monkey_ID, CELL_TYPE = 'Grans', SKIP_MONTH2 = True):
    data_expr = np.loadtxt(data_folderpath + 'data_expr.txt')
    return data_expr

def get_monkey_parameters(monkey_ID, SKIP_MONTH2 = True):
    grans_density = 6.3 * 10**8 # Ref: Chen2009reference

    #output: Mss_EGFP, sample_ratio, months
    if monkey_ID == 2:
        Mss = 7.3 * grans_density

        grans_EGFP = [5.8,5.5,5.8,4.9,6.8,5.7,5.7,6.4,9.3]
        S_EGFPs = [60900,57750,87000,36750,71400,42750,51300,48000,139500] # Ref: mmc3

        S_EGFP = np.mean(S_EGFPs[1:])
        EGFP_ratio = np.mean(grans_EGFP[1:]) * 475/679 / 100. # Ref?

        Mss_EGFP = Mss*EGFP_ratio
        sample_ratio = S_EGFP/Mss_EGFP

        months = [3, 9, 15, 29, 40, 52, 64, 81, 103]

    elif monkey_ID == 3:
        Mss = 5. * grans_density

        grans_EGFP = [3.8,6.4,6.1,7.2,6.8,7.2,6.1,6.8,7.]
        months = [2, 8, 19, 25, 32, 43, 50, 56, 67]
        S_VISs_1 = [57000.,96000.,91500.,108000.,102000.,108000.,91500.,102000.,105000.]
        S_VISs_2 = [8471., 10236., 7502., 9789., 6052., 9817., 9471., 10375., 44415.]
        num_VIS_clones = [357.,233.,195.,231.,247.,198.,255.,194.,396.]
        num_QVI_clones = [290.,184.,146.,186.,193.,152.,189.,155.,286.]

        if SKIP_MONTH2:
            grans_EGFP,months,S_VISs_1,S_VISs_2,num_VIS_clones,num_QVI_clones=\
            grans_EGFP[1:],months[1:],S_VISs_1[1:],S_VISs_2[1:],num_VIS_clones[1:],num_QVI_clones[1:]
        grans_EGFP,months,S_VISs_1,S_VISs_2,num_VIS_clones,num_QVI_clones=\
            np.array(grans_EGFP),np.array(months),np.array(S_VISs_1),np.array(S_VISs_2),\
            np.array(num_VIS_clones),np.array(num_QVI_clones)
        grans_EGFP /= 100.


        QVI_in_VIS_ratios = num_QVI_clones/num_VIS_clones


        sampled_blood_QVIs = S_VISs_1*QVI_in_VIS_ratios

         # about 77000 QVI cells in step 1

        S_VIS_avg_2 = np.mean(S_VISs_2)

        QVI_in_VIS_avg_ratio = QVI_in_VIS_ratios.mean()

        QVI_ratios_in_blood = grans_EGFP * QVI_in_VIS_ratios
        QVI_avg_ratio_in_blood = QVI_ratios_in_blood.mean()


        # From now on, EGFP_ratio is defined this way
        Mss_QVIs = Mss*QVI_ratios_in_blood
        Mss_QVI = Mss_QVIs.mean()


        S_QVIs = S_VISs_2 * QVI_in_VIS_ratios
        S_QVI_avg = S_QVIs.mean()

        #print (S_QVIs/Mss_QVIs).mean(), np.median(S_QVIs/Mss_QVIs)
        #quit()


        #print np.median(S_VISs_2), np.median(S_QVIs)

        sample_ratios = S_QVIs/Mss_QVI

        sample_ratio_avg = sample_ratios.mean()


        Mss_EGFP, sample_ratio = Mss_QVI, sample_ratio_avg


        #print QVI_in_VIS_avg_ratio

        #print 'Mss_QVI:', Mss_QVI

        #print sample_ratios, sample_ratio_avg
        #print QVI_ratios_in_blood, QVI_avg_ratio_in_blood
        #print S_QVI_avg, QVI_avg_ratio_in_blood, sample_ratio_avg

    elif monkey_ID == 4:
        Mss = 6. * grans_density

        grans_EGFP = [11, 11.3, 10.5, 10.4, 9.2, 9.7]
        months = [2, 8, 13, 19, 32, 38]
        S_EGFPs = [165000,169500,157500,156000,138000,145500]

        S_EGFP = np.mean(S_EGFPs[1:])
        EGFP_ratio = np.mean(grans_EGFP[1:]) * 1809/2213 / 100.

        Mss_EGFP = Mss*EGFP_ratio
        sample_ratio = S_EGFP/Mss_EGFP

    else:
        print 'No such monkey!'
        return None
    #print EGFP_ratio
    return [Mss_EGFP, sample_ratio, months]

def my_euclid(v1, v2, CONST_PENALTY = 1.):

    v1 = np.array(v1)
    v2 = np.array(v2)

    masks_val = (v1 == v1) & (v2 == v2)
    if sum(masks_val) == 0:
        return CONST_PENALTY * len(v1)

    v1_val = v1[masks_val]
    v2_val = v2[masks_val]

    penalty = (len(v1) - sum(masks_val)) * CONST_PENALTY

    return euclidean(v1_val, v2_val) + penalty

def my_euclid_(v1, v2, CONST_PENALTY = 1.):

    v1 = np.array(v1)
    v2 = np.array(v2)

    masks_val1 = (v1 == v1)
    masks_val2 = (v2 == v2)

    v1[masks_val1==0] = 0
    v2[masks_val2==0] = 0

    print v1, v2

    return euclidean(v1, v2)

## TODO: less sure
# Transition rate

# 0.185, ~0.88, ~1.28, ~2.5 ??

#Mss_EGFP = 5000.


def get_A(L, mu3, Mss_EGFP):
    # input: 1-by-1, 1-by-1
    # output: 1-by-1
    return Mss_EGFP * mu3 / (2**L)

def get_L(Ass_EGFP, mu3, Mss_EGFP):
    return np.log2(Mss_EGFP * mu3/Ass_EGFP)

def get_gi_s(C,H_EGFP):
    lamb = 1 - C/float(H_EGFP)

    ni_s = np.random.geometric(p=1-lamb, size=C)
    fi_s = ni_s/float(sum(ni_s))
    return fi_s

def get_sample_days_inds(sample_months, dt):
    return [int(x) for x in (np.array(sample_months)*30/dt)]

def get_HSC_diffs(A, gi_s, Tt):
    # input: 1-by-1, n-by-1
    # output: e-by-2 tuple (clone ind, diff time)
    # e is the number of events

    curr_t = 0
    events = []
    num_clones = len(gi_s)

    while curr_t < Tt:

        rate_tot = A
        curr_t += np.log(1/(1-rd.random()))/rate_tot # random() can reach 0.0

        test_num = rd.random()

        for i in range(num_clones):
            if test_num < gi_s[i]:
                events.append((i, curr_t))
                break
            else:
                test_num -= gi_s[i]

    return events[:-1] #give up the last one, which exceeds Tt

def get_m_burst(r2, L, mu3, omg, dt, Tb):
    # input: 1-by-1, 1-by-1, 1-by-1, 1-by-1, 1-by-1
    # output: t-by-1
    # t = Tb/dt

    ts_p = np.linspace(0,Tb,num=np.round(Tb/dt))

    L_thres = 15
    if L - 1 < L_thres:
        tmp = np.log(math.factorial(L-1))
    else:
        n = L-1
        tmp = n*np.log(n) - n + np.log(n*(1+4*n*(1+2*n)))/6. + np.log(np.pi)/2.
    #print L-1, np.log(math.factorial(L-1)), tmp

    nLmin1_log = np.log(2*r2*ts_p) * (L-1) - tmp - r2 * ts_p
    nLmin1 = np.exp(nLmin1_log)


    # TODO: I'm a genius! test this ? number
    # Use log and then exp to avoid computational troubles
    # For L-1 < ?, use directly
    # For L-1 >= ?, use Ramanujan's formula

    nL_one = np.convolve(nLmin1, np.exp(-ts_p*omg))*2*r2*dt # for r < omg
    #nL_one_tmp = [(2*r2/(r2-omg)) ** L * gammainc(L, (r2-omg)*t)\
    #                * np.exp(-omg*t) for t in ts_p]

    #import matplotlib.pyplot as plt
    #plt.plot(nL_one,'g')
    #plt.plot(nL_one_tmp,'b')
    #print len(nL_one), len(nL_one_tmp)
    #plt.show()

    #print np.exp(-mu3 * ts_p)
    ms_conv = np.convolve(nL_one, np.exp(-mu3 * ts_p))*omg*dt
    m_burst = ms_conv[:len(ts_p)]

    return m_burst

#m_burst = get_m_burst(r2=10., L=5, mu3=1., omg=5., dt=0.1, Tb=50)

def get_mi_s(events, m_burst, clone_num, Tt, dt):
    # input: e-by-2, t-by-1
    # output: n-by-t
    t_num = int(Tt/dt)
    mi_s = np.zeros((clone_num, t_num))

    for clone_ind, inj_time in events:
        ind_t = int(inj_time/dt)
        curr_len = min(t_num-ind_t, len(m_burst))

        mi_s[clone_ind, ind_t:ind_t+curr_len] += m_burst[:curr_len]
    return mi_s

def get_samples(mi_s, sample_days_inds, sample_ratio):
    # input: n-by-t, s-by-1, 1-by-1
    # output: n-by-s
    num_clones, num_t = mi_s.shape
    num_samples = len(sample_days_inds)
    samples = np.zeros((num_clones, num_samples))

    tmp = []
    for j in range(num_samples):
        M_EGFP = mi_s[:,sample_days_inds[j]].sum()
        # tmp.append(M_EGFP)
        # continue

        S_EGFP_theo = int(M_EGFP*sample_ratio)
        #print S_EGFP_theo

        for i in range(num_clones):
            if S_EGFP_theo <= 0:
                samples[i,j] = 0
            else:
                samples[i,j] = binomial(S_EGFP_theo,
                        mi_s[i,sample_days_inds[j]]/float(M_EGFP))
        S_EGFP_simu = sum(samples[:,j])
        if S_EGFP_simu != 0:
            samples[:,j] /= float(S_EGFP_simu)
    # print np.mean(tmp)
    return samples

def simulate_freq(r2, L, mu3, omg, gi_s, dt, Tt, Tb,
                  Mss_EGFP, sample_days_inds, sample_ratio):

    # 1.Monte Carlo on HSC differentiation events
    A = get_A(
            #r2=r2,
            L=L,
            mu3=mu3,
            Mss_EGFP=Mss_EGFP)
    HSC_diffs = get_HSC_diffs(
            A=A,
            gi_s=gi_s,
            Tt=Tt)

    clone_num = len(gi_s)

    # 2.Pre-calculate burst on r2, mu3, and dt, Tb
    m_burst = get_m_burst(r2=r2,
                          L=L,
                          mu3=mu3,
                          omg=omg,
                          dt=dt,
                          Tb=Tb)
    mi_s = get_mi_s(events=HSC_diffs,
                    m_burst=m_burst,
                    clone_num=clone_num,
                    Tt=Tt,
                    dt=dt)

    # 3.sample at times
    samples = get_samples(mi_s=mi_s,
                          sample_days_inds=sample_days_inds,
                          sample_ratio=sample_ratio)

    return samples

def get_freq(samples):
    row, col = samples.shape
    sums = np.zeros((col,))
    for j in range(col):
        sums[j] = samples[:,j].sum()
    samples_freq = np.divide(samples, sums)
    return samples_freq

def get_cumu(samples):
    row, col = samples.shape
    samples_cumu = samples
    for i in range(1,row):
        samples_cumu[i,:] += samples_cumu[i-1,:]
    return samples_cumu

def get_fz_stats(data_mat, get_clone_cnt=False):
    ################
    # Input: sample matrix
    # Output: for each z_i, get the avg & std abundances
    ################
    row, col = data_mat.shape

    zi, avgs = np.zeros((row,)), np.zeros((row,))
    for i in range(row):
        zi[i] = sum(data_mat[i,:]==0)
        avgs[i] = data_mat[i,:].mean()

    avgj_avg, avgj_std = np.zeros((col,)), np.zeros((col,))

    avgj_cnt = np.zeros((col,))
    for j in range(col):# zi=0,1,2,...,J-1
        avgj_avg[j] = avgs[zi==j].mean()
        avgj_std[j] = avgs[zi==j].std()

        avgj_cnt[j] = len(avgs[zi==j])

    if not get_clone_cnt:
        return [avgj_avg, avgj_std]
    else:
        return [avgj_avg, avgj_std, avgj_cnt]

def get_fi_all(data_mat):
    row, col = data_mat.shape

    fi_s = [[] for _ in range(col)]
    for i in range(row):
        zi = sum(data_mat[i,:]==0)
        avg = data_mat[i,:].mean()

        if zi == col:
            continue

        fi_s[int(zi)].append(avg)
    fi_s = np.array(fi_s)
    return fi_s

from scipy.stats import ttest_ind
def two_sample_test(data_expr, data_theo):
    row1, col1 = data_expr.shape
    row2, col2 = data_theo.shape
    if not (col1==col2):
        return None
    fi_s_expr = get_fi_all(data_expr)
    fi_s_theo = get_fi_all(data_theo)

    all_p_vals = []
    for j in range(1,col1): # do not include z=0
        curr_p = ttest_ind(fi_s_expr[j], fi_s_theo[j])[1]
        all_p_vals.append(curr_p)
    return all_p_vals

def resample(samples_orig, resample_num):
    row, col = samples_orig.shape
    inds_resample = rd.sample(range(col), resample_num)
    return samples_orig[:,inds_resample]

def get_clone_num_in_each_bin(data_mat, bins=None):
    row, col = data_mat.shape

    avgs = []
    for i in range(row):
        curr_avg = data_mat[i,:].mean()
        avgs.append(curr_avg)
    avgs = np.log(avgs)

    if bins==None:
        num_bins = 7
        bins = np.linspace(np.log(0.0002), np.log(0.0122), num=num_bins+1)
    avgs_in_each_bin = [[] for _ in range(num_bins)]
    for i in range(len(avgs)):
        for j in range(num_bins):
            if avgs[i] >= bins[j] and avgs[i] <= bins[j+1]:
                avgs_in_each_bin[j].append(avgs[i])

    all_num_clones = []
    for j in range(num_bins):
        all_num_clones.append(len(avgs_in_each_bin[j]))
    return all_num_clones

def get_avgstd_stats(data_mat, num_bins=7, use_linear_gap = True,
                     use_log_scale = True, use_intermediate_sizes = True):
    row, col = data_mat.shape

    avgs, stds, zi_s = [], [], []
    for i in range(row):
        avgs.append(data_mat[i,:].mean())
        stds.append(data_mat[i,:].std())
        zi_s.append(sum(data_mat[i,:]==0))

    zi_s = np.array(zi_s)
    avgs = np.log(avgs)

    if use_log_scale:
        stds = np.log(stds)
    else:
        stds = np.array(stds)

    if use_intermediate_sizes:
        avgs = avgs[zi_s > 0]
        stds = stds[zi_s > 0]
    row = len(avgs)


    if use_linear_gap:
        bins = np.linspace(np.log(0.0002), np.max(avgs), num=num_bins+1)
    else:
        avgs_sorted = np.sort(avgs)

        inds_gap = int(round(row/num_bins))

        inds = [0]
        for k in range(num_bins):
            inds.append(inds[-1] + inds_gap)

        inds[-1] = min(inds[-1], row-1)

        bins = avgs_sorted[inds]

    avgs_in_each_bin = [[] for _ in range(num_bins)]
    stds_in_each_bin = [[] for _ in range(num_bins)]
    for i in range(row):
        for j in range(num_bins):
            if avgs[i] >= bins[j] and avgs[i] <= bins[j+1]:
                avgs_in_each_bin[j].append(avgs[i])
                stds_in_each_bin[j].append(stds[i])

    avgs_avg_in_each_bin = []
    stds_avg_in_each_bin = []
    for j in range(num_bins):
        avgs_avg_in_each_bin.append(np.mean(avgs_in_each_bin[j]))
        stds_avg_in_each_bin.append(np.mean(stds_in_each_bin[j]))
    return [avgs_avg_in_each_bin, stds_avg_in_each_bin]


def get_Ri_stats(data_mat, t=1, num_bins=7, use_log_scale = True,
                 use_intermediate_sizes = True):
    row, col = data_mat.shape
    avgs = []
    Ri_s = []
    for i in range(row):
        x = data_mat[i,:]
        avgs.append(np.mean(x))

        curr_Ri = np.corrcoef(x[0:len(x)-t], x[t:len(x)])[0,1]
        Ri_s.append(curr_Ri)
    avgs = np.array(avgs)
    Ri_s = np.array(Ri_s)

    avgs = avgs[Ri_s==Ri_s]
    Ri_s = Ri_s[Ri_s==Ri_s]

    if use_log_scale:
        avgs = np.log(avgs)

    if use_intermediate_sizes:
        bins = np.linspace(np.log(0.0002), np.log(0.0122), num=num_bins+1)
    else:
        bins = np.linspace(np.log(0.0002), np.max(avgs), num=num_bins+1)


    avgs_in_each_bin = [[] for _ in range(num_bins)]
    Ri_s_in_each_bin = [[] for _ in range(num_bins)]
    for i in range(len(avgs)):
        for j in range(num_bins):
            if avgs[i] >= bins[j] and avgs[i] <= bins[j+1]:
                avgs_in_each_bin[j].append(avgs[i])
                Ri_s_in_each_bin[j].append(Ri_s[i])

    avgs_avg_in_each_bin = []
    Ri_s_avg_in_each_bin = []
    for j in range(num_bins):
        avgs_avg_in_each_bin.append(np.mean(avgs_in_each_bin[j]))
        Ri_s_avg_in_each_bin.append(np.mean(Ri_s_in_each_bin[j]))

    return [avgs_avg_in_each_bin, Ri_s_avg_in_each_bin]