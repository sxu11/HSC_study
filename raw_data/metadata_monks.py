
#import params_general as pg
import numpy as np

def get_meta_data(monkey_ID):
    grans_density = 6.3 * 10**8

    if monkey_ID == 3:
        Mss = 5. * grans_density
        months = np.array([2, 8, 19, 25, 32, 43, 50, 56, 67])
        grans_EGFP_ratio = np.array([3.8,6.4,6.1,7.2,6.8,7.2,6.1,6.8,7.])/100.

        S_VISs_1 = np.array([57000.,96000.,91500.,108000.,102000.,108000.,91500.,102000.,105000.])
        S_VISs_2 = np.array([8471., 10236., 7502., 9789., 6052., 9817., 9471., 10375., 44415.])

        num_VIS_clones = np.array([357.,233.,195.,231.,247.,198.,255.,194.,396.])
        num_QVI_clones = np.array([290.,184.,146.,186.,193.,152.,189.,155.,286.])

        # EGFP+/total is directly measurable
        # VIS/QVI is approximated by clone-wise VIS/QVI

        DNAs_grans = [10.] * len(months)

    elif monkey_ID == 2:
        Mss = 7.3 * grans_density
        months = np.array([3, 9, 15, 29, 40, 52, 64, 81, 103])
        grans_EGFP_ratio = np.array([5.8,5.5,5.8,4.9,6.8,5.7,5.7,6.4,9.3])/100.
        S_VISs_1 = np.array([60900.,57750.,87000.,36750.,71400.,42750.,51300.,48000.,139500.])
        S_VISs_2 = np.array([4732., 4072., 4309., 2949., 3794., 3362., 11170., 5430., 7777.])
        num_VIS_clones = np.array([157., 193., 268., 176., 213., 169., 183., 151., 182.])
        num_QVI_clones = np.array([128., 155., 217., 126., 147., 114., 141., 129., 156.])

        # EGFP+/total is directly measurable
        # VIS/QVI is approximated by clone-wise VIS/QVI

        DNAs_grans = [7., 7., 10., 5., 7., 5., 6., 5., 10.] # checked

    elif monkey_ID == 4:
        Mss = 6. * grans_density
        months = np.array([2, 8, 13, 19, 32, 38])

        grans_EGFP_ratio = np.array([11, 11.3, 10.5, 10.4, 9.2, 9.7])/100.
        S_VISs_1 = np.array([165000.,169500.,157500.,156000.,138000.,145500.])
        S_VISs_2 = np.array([8531., 12860., 10451., 9935., 8214., 9793.])
        num_VIS_clones = np.array([505.,614.,689.,392.,523.,942.])
        num_QVI_clones = np.array([435.,516.,570.,324.,449.,791.])


        # EGFP+/total is directly measurable
        # VIS/QVI is approximated by clone-wise VIS/QVI

        DNAs_grans = [10.] * len(months)


    grans = [DNAs_grans[i]*150000 for i in range(len(months))]


    grans_EGFP = [grans[i]*grans_EGFP_ratio[i] for i in range(len(months))]

    QVI_in_VIS_ratios = num_QVI_clones/num_VIS_clones
    S_QVIs = S_VISs_2 * QVI_in_VIS_ratios

    print np.mean(grans_EGFP_ratio)* QVI_in_VIS_ratios.mean()

    return grans, grans_EGFP, S_QVIs, months

get_meta_data(4)
