

import pandas as pd
pd.set_option('display.width', 3000)
pd.set_option('display.max_columns', 3000)

df = pd.read_csv('abunds/zh33.csv')

cell_types = ['T', 'B', 'Mon', 'Gr', 'NK']
times = ['1m', '2m', '3m', '4.5m', '6.5m', '9.5m']
ts = [float(x[:-1]) for x in times]

cols = {}
for time in times:
    cols[time] = []
    for col in df.columns.values:
        if (time in col) and ('node' not in col):
            cols[time].append(col)

cols['6.5m'] = ['6.5m T', '6.5m B', '6.5m NK overall', '6.5m Mon', '6.5m Gr']

def alph_diversity(arr):
    return sum(arr>0)
def gamm_diversity(arr):
    return 5
def beta_diversity(arr):
    return gamm_diversity(arr) - alph_diversity(arr)

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

'''https://ecopy.readthedocs.io/en/latest/diversity.html'''
import ecopy
# for time in times:
#     print cols[time]
D_alpha, D_beta, D_gamma = ecopy.div_partition(df[cols['6.5m']].transpose(),
                                               method='spRich')
# print D_alpha, D_beta, D_gamma
#
#
# df_tmp = df[cols['6.5m']]
# df_tmp = df_tmp[(df_tmp.T != 0).any()]
# print df_tmp.shape
# quit()
#
# tmps = []
# for i in range(5):
#     tmps.append(sum(df_tmp[cols['6.5m'][i]]>0))
# print sum(tmps)/5.
# quit()

alphs, betas, gammas = [], [], []
for time in times:
    D_alpha, D_beta, D_gamma = ecopy.div_partition(df[cols[time]].transpose(),
                                                   method='spRich')
    alphs.append(D_alpha)
    betas.append(D_beta)
    gammas.append(D_gamma)

matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)

fig, ax1 = plt.subplots()
ax1.plot(ts, alphs, label=r'$\alpha$')
ax1.plot(ts, gammas, label=r'$\gamma$')
ax1.set_xlabel('month', fontsize=16)
ax1.set_ylabel(r'$\alpha$ and $\gamma$ diversity', fontsize=16)

plt.legend(loc=(.8,.8), fontsize=14)
#


ax2 = ax1.twinx()

ax2.plot(ts, betas, 'k--', label=r'$\beta$')
ax2.set_ylabel(r'$\beta$ diversity', fontsize=16)
plt.legend(loc=(.8,.72), fontsize=14)
plt.tight_layout()
plt.show()

# T_cols, B_cols, Mon_cols, Gr_cols, NK_cols = [], [], [], [], []
# for col in df.columns.values:
#     if 'T' in col:
#         T_cols.append(col)
#     elif 'B' in col:
#         B_cols.append(col)
#     elif 'Mon' in col:
#         Mon_cols.append(col)
#     elif 'NK' in col:
#         NK_cols.append(col)
#     elif 'Gr' in col:
#         Gr_cols.append(col)
#
# print T_cols, B_cols, Mon_cols, Gr_cols, NK_cols
# df_T = df[T_cols[:-1]]
# df_B = df[B_cols[:-1]]
#
# print df[Mon_cols]