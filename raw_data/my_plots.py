
import matplotlib.pyplot as plt
import matplotlib as mpl

def use_sci_nota(usesci_x=False, usesci_y=False):

    ax = plt.gca()

    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
    if usesci_x:
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
    if usesci_y:
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


# [FONTSIZE_MATH,FONTSIZE_TEXT,FONTSIZE_TICK,FONTSIZE_LEGD] = [36,32,28,27] #[36,32,28,24] #[50,45,35,35]
# [MARGIN_LEFT,MARGIN_RIGHT,MARGIN_BOTTOM,MARGIN_TOP] = [.18,.9,.15,.9]#[.27,.93,.15,.9]#[.21,.93,.15,.9]# #

[FONTSIZE_MATH,FONTSIZE_TEXT,FONTSIZE_TICK,FONTSIZE_LEGD] = [30,30,24,20] # 30,30,24,24
[MARGIN_LEFT,MARGIN_RIGHT,MARGIN_BOTTOM,MARGIN_TOP] = [.18,.87,.15,.9]

def config_plot(xlab='', ylab='', legend_loc=(.8,.7), legend_prop=18):


    plt.rc('font', size=FONTSIZE_TICK, weight='bold')
    plt.tick_params(labelsize=FONTSIZE_TICK, width=1.5, which='both', direction='out')
    plt.rc('text', usetex=True)
    mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

    x_is_math = '$' in xlab
    y_is_math = '$' in ylab

    #ax = plt.gca()
    #ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))

    plt.legend(bbox_to_anchor=(legend_loc),bbox_transform=plt.gcf().transFigure,
               fontsize=FONTSIZE_LEGD, prop={'size':FONTSIZE_LEGD},frameon=False)


    plt.xlabel(xlab, fontsize=FONTSIZE_MATH if x_is_math else FONTSIZE_TEXT)
    plt.ylabel(ylab, fontsize=FONTSIZE_MATH if y_is_math else FONTSIZE_TEXT)
    plt.subplots_adjust(left=MARGIN_LEFT, bottom=MARGIN_BOTTOM,
                        right=MARGIN_RIGHT, top=MARGIN_TOP)