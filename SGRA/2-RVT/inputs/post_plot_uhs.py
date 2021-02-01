# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 23:44:11 2021

@author: FELi
plot surface UHS, in comparison with bedrock UHS & surface spectrum using the multiplication appraoch
"""
import numpy as np
from math import exp, log
from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns


## rp of interest
list_rp = [2475, 5000, 10000]

fit_folder = r'C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\20368492 National Grid\SGRA\2-RVT/'
fname = fit_folder+'uhs_surface.csv' #'CMS_NGA_east.csv'
uhs_all = np.genfromtxt(fname, delimiter=',', dtype=None, names=True)

### plot UHS
uhs_all = np.array(uhs_all.tolist())
# print(type(uhs_all))
# print(uhs_all)
# print('uhs_all[:,0]=', uhs_all[:,0])
# N = len(list_rp)
# temp = sns.color_palette("hls", N)
# list_colors = temp+ temp +temp

list_ltypes = ["solid",#"solid","solid",
               "dashed",#"dashed","dashed",
               "dotted",]#"dotted","dotted",] 

# print("list_colors", list_colors)

uhs_folder = r"C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\Bedrock PSHA results\Salem\UHRS/"
## surface spectra using traditional approach

traditional_file = fit_folder+'RSP_surface_traditional.csv'
trad_data = np.genfromtxt(traditional_file, delimiter=',', dtype=None, names=True)


# trad_data = np.array(trad_data.tolist())

for ii in range(len(list_rp)):
    irp = list_rp[ii]
    
    fig, ax0 = plt.subplots(nrows=1, ncols=1, sharex=True,
                                                figsize=(7, 6))
    icount = 0
    
    uhs_file = uhs_folder+'UHS_'+str(irp)+'yr.txt'
    uhs_data = np.genfromtxt(uhs_file, delimiter=' ', dtype=None)
    # print(uhs_data)
    # print(type(uhs_data))
    # UHS_periods = uhs_data[:,0]
    # UHS_sa = uhs_data[:,1]
    ax0.semilogx(uhs_data[:,0],uhs_data[:,1],
                 color = 'k', #list_colors[icount], 
                 # marker='o', 
                 linestyle=list_ltypes[icount],
                 linewidth=1, #markersize=12,
                 label= f'bedrock-{list_rp[icount]}yr')
        
    
    icount+=1
    ax0.semilogx(uhs_all[:,0],uhs_all[:,ii+1],
                 color = 'k', # list_colors[icount], 
                 # marker='o', 
                 linestyle=list_ltypes[icount],
                 linewidth=1, #markersize=12,
                 label= f'surface-{list_rp[ii]}yr')
        
    ### deterministic AF (as the envelop of all CMS)
    icount+=1

    ax0.semilogx(trad_data["period_s"],trad_data[f"yr{irp}"],
                 color = 'k', # list_colors[icount], 
                 # marker='o', 
                 linestyle=list_ltypes[icount],
                 linewidth=1, #markersize=12,
                 label= f'surface-{list_rp[ii]}yr-x')
    
    ax0.legend() 
    
    ####----
    # ax0.set_ylim([0.0001,0.1])
    ax0.set(xlabel='Period (s)', ylabel=f'SA (g)',
                    title=f'{irp}-year')
    ax0.grid()
    
    fig.savefig(fit_folder+f'Salem_uhs_surface_{irp}yr.png', dpi=300)
    plt.show() 

