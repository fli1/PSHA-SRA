# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 12:01:28 2021

@author: FELi
plot a selected param, e.g. csr, over bedrock motions

"""

import numpy as np
import matplotlib.pyplot as plt
# from math import median
from statistics import median
# from math import exp
import matplotlib.backends.backend_pdf



def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


rvt_output_folder = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\SGRA\2-RVT/'

cms_initial = 'Salem_cms_initial.csv'
# param = "Amplification Factor"
# param_file = 'RVTmulti-0.00 ft (Outcrop (2A)) from Bedrock (Outcrop (2A))-specRatio.csv'
# coeff = 1.0
param = "CSR"
param_file = 'RVT1-profile-stressRatio.csv'
coeff = 0.65

param_data = np.genfromtxt(rvt_output_folder+param_file, delimiter=',', dtype=None,#names=True, 
                        skip_header=3)
cmsdata  = np.genfromtxt(rvt_output_folder+'inputs/'+cms_initial, delimiter=',', dtype=None, names=True)

# print(len(cmsdata))
cmsdata = np.array(cmsdata.tolist())
cms_periods = cmsdata[:,0]
no_periods = len(cms_periods)
no_cms = len(cmsdata[0,:])-1
idx = np.where(cms_periods==0.01)
pga_bedrock = cmsdata[idx,1:].flatten()
print(pga_bedrock)

# print(type(param_data))
# print(type(param_data[0]))
# print(param_data[0])

param_data = np.array(param_data.tolist())
# param_periods = param_data[:,0]

no_params = len(param_data[0,:])-3
no_sims = no_params/no_cms
no_depths = len(param_data[:,0])

ii = 0

pdf = matplotlib.backends.backend_pdf.PdfPages(rvt_output_folder+f"{param}_output.pdf")

for ii in range(no_depths):
    idepth = param_data[ii,0]

    print(idepth)
    
    # idx,param_p = find_nearest(param_periods, ip)    
    fig, ax0 = plt.subplots(nrows=1, ncols=1, sharex=True,
                                                figsize=(7, 6))
    median_out = []
    sa_bedrock_all = []
    param_all = []
    
    for jj in range(no_cms):
        temp_idx = range(jj+1,no_params+1, no_cms)
        param_datasub = param_data[ii,temp_idx]*coeff
        
        median_out.append(median(param_datasub))
        
        temp_x = np.repeat(pga_bedrock[jj],no_sims)
        
        # print(temp_x)
        # print(param_datasub)

        ax0.semilogx(temp_x, param_datasub, "ko")
        sa_bedrock_all.append(temp_x)
        param_all.append(param_datasub)
    
    
    sa_bedrock_all = np.asarray(sa_bedrock_all).flatten()
    param_all = np.asarray(param_all).flatten()
    # print(type(sa_bedrock_all))
    # print(type(param_all))
    # print(param_all)
    p = sa_bedrock_all.argsort()

    x=  sa_bedrock_all[p]
    y = np.asarray(param_all)[p]
    # ax0.semilogx(x,y, "k-")

    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    xp = np.logspace(-4, 0, 50)
    ax0.semilogx(xp,p(xp), "r-")
    
    # ax0.set_ylim([0,4])
    ax0.set(xlabel='PGA_bedrock (g)', ylabel=f'{param}',
                    title=f'Depth={idepth}ft')
    ax0.grid()
    
    # fig.savefig(rvt_output_folder+f'{param}_Salem_{ii}_{idepth}ft.png', dpi=300)
    plt.show()
    pdf.savefig( fig )
    
pdf.close()              