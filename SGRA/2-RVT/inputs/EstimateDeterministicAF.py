# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 23:55:05 2021

@author: FELi
estimate the median AF for each return period
"""

import numpy as np
# ### plot UHS
# uhs_folder = r"C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\Bedrock PSHA results\Salem\UHRS/"
# list_rp = [10, 50, 100, 250, 475, 975, 2475, 5000, 10000]
# # for irp in list_rp:
# irp=50
# uhs_file = uhs_folder+'UHS_'+str(irp)+'yr.txt'
# uhs_data = np.genfromtxt(uhs_file, delimiter=' ', dtype=None)
# print(uhs_data)
# # print(type(uhs_data))
# UHS_periods = uhs_data[:,0]
# UHS_sa = uhs_data[:,1]
# print(UHS_periods)

fit_folder = r'C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\20368492 National Grid\SGRA\2-RVT/'
filename = 'RVT1-0.00 ft (Outcrop (2A))-respSpec.csv'

ifile = fit_folder+filename
param_data = np.genfromtxt(ifile, delimiter=',', dtype=None, skip_header=3)#,names=True)
# print(param_data)

list_rp = [50, 100, 250, 475, 975, 2475, 5000, 10000]

list_target = ['uhs', 'uhs', 'cms','cms','cms',
                    'cms','cms','cms']
list_periods = [[0.2], [0.2], [0.2,1], [0.2,1,5], [0.2,1,5], 
                [0.2,1,5], [0.2,1,5], [0.2,1,5], ]

no_cms = 19
no_params = no_cms*100
temp_idx = 0
case = 'period_s'
dataout = param_data[:,0]#.tolist()
for ii in range(len(list_rp)):
    irp = list_rp[ii]
    itarget = list_target[ii]
    ip = list_periods[ii]
    print(ip)
    print(len(ip))
    if itarget=='uhs':
        # temp_idx = 
        jj_cols = range(temp_idx+1,no_params+1, no_cms)
        param_datasub = param_data[:,jj_cols]
        AF_median = np.median(param_datasub, axis=1)
        case += f',yr{irp}'
        temp_idx += 1
    if itarget=='cms':
        # temp_idx = 
        jj_cols = range(temp_idx+1,no_params+1, no_cms)
        param_datasub = param_data[:,jj_cols]
        AF_median = np.median(param_datasub, axis=1)
        case += f',yr{irp}'
        temp_idx += 1
        
        # AF_median = 
        for iip in range(1,len(ip)):
            jj_cols = range(temp_idx+1,no_params+1, no_cms)
            param_datasub = param_data[:,jj_cols]
            temp = np.median(param_datasub, axis=1)
            AF_median = np.maximum(AF_median, temp)
            # case += f'yr{irp}_{itarget}'
            temp_idx += 1
    # print(AF_median)
    
    dataout = np.vstack((dataout,AF_median))
    # dataout.append(AF_median.tolist())

print(dataout)
dataout = np.transpose(dataout)
fname = fit_folder+'RSP_surface_traditional.csv' #'CMS_NGA_east.csv'
with open(fname, 'w') as f:    
    # np.savetxt(f, case, delimiter=',')
    f.write(case + '\n')
with open(fname, 'ba') as f:
    np.savetxt(f, dataout, delimiter=',')       
        
