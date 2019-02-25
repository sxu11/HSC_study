
# For Monkey 3
months =  [2, 8, 19, 25, 32, 43, 50, 56, 67]
DNAs = [10.] * len(months)
monkey_weight = 5.
grans_EGFP = [3.8,6.4,6.1,7.2,6.8,7.2,6.1,6.8,7.]
pbmc_EGFP = [1.9, 4.8, 6.2, 5.8, 5.8, 6.4, 6.1, 6.5, 6.6]
EGFP_ratio = 0.05 # sum(grans_EGFP[1:])/len(grans_EGFP[1:])/100.

thres_large = 0.0122
thres_small = 0.0005