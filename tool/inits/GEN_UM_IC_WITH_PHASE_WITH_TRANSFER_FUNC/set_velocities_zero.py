#!/usr/bin/env python3.7

##################################################################################################
# This script is used for setting the velocity-dependent transfer fucntion (TF) to 1e-8 for TF   #
# file generated by CAMB(planck_2018_transfer_out.dat) and axsionCAMB(planck_2018_axion_transfer #
# _out.dat). Be sure to prepare the input transfer function files and specify their name!! The   #
# output file planck_2018_transfer_out_no_vel.dat/planck_2018_axion_transfer_out.dat will be     #
# produced accordingly(corresponding to CAMB/axionCAMB input).                                   #
##################################################################################################

import sys
import numpy as np
import pandas as pd

if len(sys.argv) != 2:
    raise RuntimeError("Input Error: Please select do axion or not (T/F)!")

if sys.argv[1] == "F":
    # no axion
    df_trans = pd.read_csv('planck_2018_transfer_out.dat', delimiter='\s+', index_col=False)
    col      = list(df_trans.columns[1:])
    col.append("None")
    df_trans.columns = col
    df_trans = df_trans.iloc[:,:-1]
    df_trans['v_CDM'] = 1e-8*np.ones(df_trans['v_CDM'].shape)
    df_trans['v_b'] = 1e-8*np.ones(df_trans['v_b'].shape)
    df_trans['v_b-v_c'] = 1e-8*np.ones(df_trans['v_b-v_c'].shape)

    df_trans.to_csv('planck_2018_transfer_out_no_vel.dat',header=None, sep='\t', float_format = "%.6e", index=None)
####################################################################################################################
elif sys.argv[1] == "T":
    # with axion
    df_trans_axion = pd.read_csv('planck_2018_axion_transfer_out.dat', delimiter='\s+', index_col=False, header=None)

    # switch the CDM column with axion column
    temp                      = df_trans_axion.iloc[:,1].copy()
    df_trans_axion.iloc[:,1]  = df_trans_axion.iloc[:,6].copy()
    df_trans_axion.iloc[:,6]  = temp.copy()

    # switch the CDM column with total matter column
    temp                      = df_trans_axion.iloc[:,6].copy()
    df_trans_axion.iloc[:,6]  = df_trans_axion.iloc[:,-1].copy()
    df_trans_axion.iloc[:,-1] = temp.copy()

    # subtitute the data with negative TF by its absolute value
    df_trans_axion = abs(df_trans_axion)

    # subtitute the data with negative TF by 1e-8
    for i in range(9,13):
        df_trans_axion[i] = 1e-8*np.ones(len(df_trans_axion.iloc[:,0]))

    df_trans_axion.iloc[:,:100].to_csv('planck_2018_axion_transfer_out_no_vel.dat', \
                                        header=None, sep='\t', float_format = "%.6e", index=None)
####################################################################################################################
else:
    raise RuntimeError("Input Error: Please select do axion or not (T/F)!")
