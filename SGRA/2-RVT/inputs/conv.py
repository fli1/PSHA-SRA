# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 23:27:41 2021

@author: FELi

CONVOLUTION TO GET SURFACE SPECTRA

"""
import numpy as np
from math import exp, log
from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns

## rp of interest
list_rp = [2475, 5000, 10000]
list_aep = [1/2475., 1/5000., 1/10000.]

list_periods = [0.01,0.02,0.03,0.05,#0.07,
                0.1,0.15,
                0.2,0.25,0.3,0.4,0.5,0.75,1,
                1.5,2,3,4,5,7.5,10]

### read AF functions
fit_folder = r'C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\20368492 National Grid\SGRA\2-RVT/'
fit_file = 'fit_aprams.csv'

fit_data = np.genfromtxt(fit_folder+fit_file, delimiter=',', dtype=None,
                        names=True)
fit_data = np.array(fit_data.tolist())
fitdata_periods = fit_data[:,0]

## read hc data
hc_folder = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\Bedrock PSHA results\Salem\hazard curves/'
uhs_all = []


hc_out = []
for ip in list_periods:        
    ### read hazard curve data
    # ip = 'PGA'
    # ip=0.2
    if ip==0.01:
        hc_file = f'curves_PGA.csv'
    else:
        hc_file = f'curves_{ip}.csv'
    
    hc_data = np.genfromtxt(hc_folder+hc_file, delimiter=',', dtype=None,
                           # names=False, 
                            skip_header=0)
    # print(hc_data[0])
    # print(type(hc_data[0]))
    sa_bedrock = hc_data[0].tolist()[3:]
    aep_bedrock = hc_data[1].tolist()[3:]
    # print('aep_bedrock=',aep_bedrock)
    # print(len(aep_bedrock))
    
    ### AF function
    idx = np.where(fitdata_periods==ip)[0][0]
    f1 = fit_data[idx,1]
    f2 = fit_data[idx,2]
    f3 = fit_data[idx,3]
    sd = fit_data[idx,4]
    # print('idx=', idx, ', sd=', sd)
    # print('f1=',fit_data[0,1])
    # print('f1=',fit_data[idx,1])
    ### convolution
    sa_surface = sa_bedrock
    aep_surface_all = []
    
    for isa_surface in sa_surface:
        aep_surface = 0.
        
        # ii=0
        aep_bedrock_pre = aep_bedrock[0]
        aep_bedrock_next = aep_bedrock[1]
        
        for ii in range(len(sa_bedrock)):
            isa_bedrock = sa_bedrock[ii]
            
            af =  isa_surface/isa_bedrock
            af_median_ln = f1+f2*log((f3+isa_bedrock)/f3)
            af_sd = sd
            p_af = 1-norm(loc = af_median_ln , scale = af_sd).cdf(log(af))
            # print('p_af=',p_af,', af_median_ln=',af_median_ln,
            #       ', af_sd=',af_sd, ', log(af)=', log(af))
            # print('af =', af)
            aep_surface += p_af*((aep_bedrock_pre-aep_bedrock_next)/2)
            # ii+=1
            aep_bedrock_pre = aep_bedrock[ii]
            # print(ii, len(aep_bedrock))
            if ii>len(aep_bedrock)-3:
                aep_bedrock_next = aep_bedrock[ii]
            else:
                aep_bedrock_next = aep_bedrock[ii+2]
            # print('aep_bedrock_pre=',aep_bedrock_pre,', aep_bedrock_next=',aep_bedrock_next)
        aep_surface_all.append(aep_surface)
        # print('aep_surface=',aep_surface)
        # print(type(aep_surface))
    
    # aep_surface_all = aep_surface_all.flatten()
    # print(aep_surface_all)
    # print(type(aep_surface_all))
    ## calculate UHS
    uhs = np.interp(np.log(list_aep), np.log(aep_surface_all)[::-1], np.log(sa_surface)[::-1])
    print('uhs=', np.exp(uhs))
    print('list_aep=', list_aep)
    print('aep_surface_all=',aep_surface_all)
    print('sa_surface=',sa_surface)
    uhs_all.append(np.append(ip,np.exp(uhs)))
    
    hc_out.append(sa_surface)
    hc_out.append(aep_surface_all)
    
    ### plot
    fig, ax0 = plt.subplots(nrows=1, ncols=1, sharex=True,
                                                figsize=(7, 6))
    
    ax0.loglog(sa_surface,aep_surface_all, "k-", label=f'{ip}s surface')
    ax0.loglog(sa_bedrock,aep_bedrock, "k--", label=f'{ip}s bedrock')
    
    # dsha_sa_surface = []
    # for ii in sa_bedrock:    
    #     af_median = exp(f1+f2*log((f3+ii)/f3))
    #     dsha_sa_surface.append(af_median*ii)
    #     print('AF=', af_median, ', sa_bedrock=', ii)
        
    # # dsha_sa_surface = dsha_sa_surface.flatten()
    # ax0.loglog(dsha_sa_surface,aep_bedrock, "b--", label=f'bedrock*AF')
    ax0.legend() 
    
    ####----
    ax0.set_ylim([0.0001,0.1])
    ax0.set(xlabel='SA_surface (g)', ylabel=f'AEP',
                    title='')
    ax0.grid()
    
    fig.savefig(fit_folder+f'Salem_hc_surface_{ip}s.png', dpi=300)
    plt.show()
        

fname = fit_folder+'uhs_surface.csv' #'CMS_NGA_east.csv'
with open(fname, 'w') as f:    
    # np.savetxt(f, case, delimiter=',')
    f.write('period_s,' + '\n')
with open(fname, 'ba') as f:
    np.savetxt(f, uhs_all, delimiter=',')

