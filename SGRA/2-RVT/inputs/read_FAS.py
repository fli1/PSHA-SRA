# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 13:18:47 2021

@author: FELi
read FAS and plot from json file to verify inputs

"""
import matplotlib.pyplot as plt
import json
import numpy as np
from math import log, log10
import seaborn as sns
        
        
N = 19
list_colors = sns.color_palette("hls", N)
list_ltypes = ["solid", "dashed","solid", "dashed","solid", "dashed",
                "solid", "dashed","solid", "dashed","solid", "dashed",
                "solid", "dashed","solid", "dashed","solid", "dashed",
                "solid", "dashed","solid", "dashed","solid", "dashed"]
rvt_folder = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\SGRA\2-RVT\inputs/'
json_file = 'Salem_RVT_empirical_inputs.json'

freq_max = 50
freq_min = 0.05
freq_size = 1024
freq_spacing =1
freq = np.logspace(log10(freq_min),log10(freq_max),freq_size)
fourierAcc_out = []
fig, ax0 = plt.subplots(nrows=1, ncols=1, sharex=True,
                                                figsize=(7, 6))

ii=-1
with open(rvt_folder+json_file, 'r') as f:
    data = json.load(f)
    motion_library = data['motionLibrary']
    motions = motion_library['motions']
    for cur_motion in motions:
        ii +=1
        fourierAcc = cur_motion['fourierAcc']
        description = cur_motion['description']
        fourierAcc_out.append(fourierAcc)
        ax0.semilogx(freq, fourierAcc, 
                     color = list_colors[ii], 
                      # marker='o', 
                      linestyle=list_ltypes[ii],
                      linewidth=1, #markersize=12,
                      label=description)

ax0.set(xlabel='Freq (hz)', ylabel=f'Fourier Amp',
                    title='FAS')
ax0.grid()
ax0.legend()
fig.savefig(f'FAS_Salem_.png', dpi=300)
plt.show()


RSPAcc_out = []
fig, ax0 = plt.subplots(nrows=1, ncols=1, sharex=True,
                                                figsize=(7, 6))

ii=-1
with open(rvt_folder+json_file, 'r') as f:
    data = json.load(f)
    motion_library = data['motionLibrary']
    motions = motion_library['motions']
    for cur_motion in motions:
        ii +=1
        period = cur_motion['targetRespSpec']['period']
        Acc = cur_motion['targetRespSpec']['sa']
        description = cur_motion['description']
        RSPAcc_out.append(Acc)
        if(ii>10):
            ax0.semilogx(period, Acc, 
                         color = list_colors[ii], 
                          # marker='o', 
                          linestyle=list_ltypes[ii],
                          linewidth=1, #markersize=12,
                          label=description)

ax0.set(xlabel='Period (s)', ylabel=f'Sa (g)',
                    title='RS')
ax0.grid()
ax0.legend()
fig.savefig(f'RS_Salem_.png', dpi=300)
plt.show()

