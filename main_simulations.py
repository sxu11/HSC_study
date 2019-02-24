
import utils_simulations as rs
import numpy as np
import os.path

curr_py_ID = 1

is_test_mode = False
# Purpose 3: fix a mu3, get the best estimate of (L, r2) pair

monkey_ID = 3
samples_expr = rs.get_expr_data(monkey_ID=monkey_ID)

[Mss_EGFP, sample_ratio, months] = rs.get_monkey_parameters(monkey_ID=monkey_ID)

mu3 = 1. # fixed

if is_test_mode:
    L_s = [24.]
    r2_s = [1.]
else:
    #L_s = np.linspace(16,18.5,num=6).tolist() + np.linspace(19,28,num=46).tolist() \
    #     + np.linspace(28.5,32,num=8).tolist()

    L_s = np.linspace(19,28,num=91).tolist()
    #np.linspace(19,21.5,num=6).tolist() + np.linspace(22,25,num=31).tolist() + \
        #np.linspace(25.5,28,num=6).tolist()
    #Ass_EGFP_s = [1,10,100]
    #print (L_s)
    r2_s = np.linspace(0.1,6.,num=60).tolist()
    #np.linspace(0.1,1.,num=10).tolist() + np.linspace(1.5,6,num=10).tolist()

    #Ls_len = len(L_s)
    #L_s = L_s[curr_py_ID*2-2:curr_py_ID*2] + \
    #      L_s[Ls_len-(curr_py_ID*2):Ls_len-(curr_py_ID*2-1)+1]

Tt, dt = 70 * 30 + 1, 1. # dt is in /day

omg = 1/6.2 # Ref: Dancey1967

res_folderpath = 'data_file/'
if not os.path.exists(res_folderpath):
    os.mkdir(res_folderpath)

for L in L_s:
    for r2 in r2_s:
        res_filepath = res_folderpath + '[r2,L]='+"[%.1f,%.1f]"%(r2,L)
        if os.path.isfile(res_filepath):
            print res_filepath + ' has been processed!'
            continue

        A = rs.get_A(L=L,
                     mu3=mu3,
                     Mss_EGFP=Mss_EGFP)
        print res_filepath, 'A:', A
        gi_s = np.mean(samples_expr, axis=1)

        Tb= max((L/r2 + 1/omg + 1/mu3)*2., 80)

        # do a simulation
        HSC_diffs = rs.get_HSC_diffs(
            A=A,
            gi_s=gi_s,
            Tt=Tt)

        clone_num = len(gi_s)

        # 2.Pre-calculate burst on r2, mu3, and dt, Tb
        m_burst = rs.get_m_burst(r2=r2,
                              L=L,
                              mu3=mu3,
                              omg=omg,
                              dt=dt,
                              Tb=Tb)

        mi_s = rs.get_mi_s(events=HSC_diffs,
                        m_burst=m_burst,
                        clone_num=clone_num,
                        Tt=Tt,
                        dt=dt)

        # save the sample matrix

        mi_s = mi_s[mi_s.sum(axis=1) > 0, :]
        np.savetxt(res_filepath, mi_s, fmt="%.0f")
