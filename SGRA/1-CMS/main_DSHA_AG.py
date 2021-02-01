# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 22:26:11 2021

@author: FELi

run NGA-East DSHA
"""
from NGA_East_20200114 import runNGAEastGMM 
##from NGA_east_SiteTerm import calF760, calFv
import numpy as np
import math
from USGS_siteamp_AG import site_amp
from USGS_sigma_AG import sigmatotal
from math import log , exp, sqrt

import matplotlib
import matplotlib.pyplot as plt
from numpy import genfromtxt
import pandas as pd

import seaborn as sns

####################################################### 
################ utility ########################### 
####################################################### 
def nga_east_calc(inPeriodList,inMagnitude,inRrup,inSigmaType,inNsigma,
                  vs30):
    
    # print("inRrup ",type(inMagnitude))
    # print("inRrup ", type(inRrup))
    DataMatrix = runNGAEastGMM(inPeriodList,inMagnitude,inRrup,inSigmaType,inNsigma)
    # DataMatrix.append([thisT,medianPSA,MedianMinusNSigma,MedianPlusNSigma])
    
    ## usgs2018 sigma model
    sigma = sigmatotal(inMagnitude, vs30,
                    inPeriodList)  
    
    ##---- usgs2018 site term
    pgaRock = DataMatrix[0][1]
    fT, sT = site_amp(inPeriodList, vs30, pgaRock)  

    # 
    PSA = []
    PSAm = []
    for ii in range(0, len(inPeriodList)):
        PSAmedian = DataMatrix[ii][1]*fT[ii]
        PSA84th = math.exp(math.log(PSAmedian)+sigma[ii])
        PSA84th_4comp = DataMatrix[ii][3]*fT[ii]
        PSA.append([inPeriodList[ii], PSAmedian, PSA84th, PSA84th_4comp])
        PSAm.append(PSAmedian)

    title = 'period, median, 84th, 84th(NGA-east_sigma)'
    with open('psa_nga_east.csv', 'w') as f:
        f.write(title + '\n')
    with open('psa_nga_east.csv', 'ba') as f:
        np.savetxt(f, PSA, delimiter=',')
    
    # print("PSAmedian= ",PSAm)
    return PSAm, sigma


def baker_jayaram_correlation(T1, T2):

# Created by Feng LI, 7/30/2019 
# Compute the correlation of epsilons for the NGA ground motion models#%
#% The function is strictly emperical, fitted over the range the range 0.01s <= T1, T2 <= 10s
#%#% Documentation is provided in the following document:
#% Baker, J.W. and Jayaram, N., "Correlation of spectral acceleration values from NGA ground 
#% motion models," Earthquake Spectra, (in review).
#
#% INPUT
#%#%   T1, T2      = The two periods of interest. The periods may be equal,
#%                 and there is no restriction on which one is larger.
#%
#% OUTPUT
#%#%   rho         = The predicted correlation coefficient

    T_min = min(T1, T2)
    T_max = max(T1, T2)
    
    C1 = (1-math.cos(math.pi/2 - math.log(T_max/max(T_min, 0.109)) * 0.366 ))
    if T_max < 0.2 :
        C2 = 1 - 0.105*(1 - 1./(1+math.exp(100*T_max-5)))*(T_max-T_min)/(T_max-0.0099)
    if T_max < 0.109:
        C3 = C2
    else:
        C3 = C1
    
    C4 = C1 + 0.5 * (math.sqrt(C3) - C3) * (1 + math.cos(math.pi*(T_min)/(0.109)))

    if T_max <= 0.109:
        rho = C2
    elif T_min > 0.109:
        rho = C1
    elif T_max < 0.2:
        rho = min(C2, C4)
    else:
        rho = C4
    
    return(rho)

def CMS_calc(list_Tstar, inPeriodList, 
            list_mag,list_R, 
            UHS_periods, UHS_sa,
            outcsvfilename):
    #% Ground motion prediction model calcs
    case_no = len(list_Tstar)*2+2
    period_no = len(inPeriodList)
    totalpoints = case_no*period_no
    cms_out = np.array(range(totalpoints)).reshape(period_no,case_no).astype(np.float32)
    # print("csm_out", cms_out)
    case = 'period, target, '
    icase = -1
    for kk in range(len(list_Tstar)):
        print(list_mag[kk])
        iM = round(float(list_mag[kk]),2)
        iR = round(float(list_R[kk]),2)
        PSAm, sigma = nga_east_calc(inPeriodList, iM, iR,
                                 inSigmaType, inNsigma, vs30)
    
        iT = list_Tstar[kk]
        #% epsilon at T_star
        sa1 = np.interp(iT, UHS_periods, UHS_sa)
        sa2 = np.interp(iT, inPeriodList, PSAm)
        ss = np.interp(iT, inPeriodList, sigma)
        
        epsilon = (log(sa1)-log(sa2))/ss; 
        
        cms = []
        cms_sigma = []
        icase +=1
        
        for ii in range(len(inPeriodList)):
            ip = inPeriodList[ii]
            
            target_sa = np.interp(ip, UHS_periods, UHS_sa)
            
            sa_median = PSAm[ii]
            ss_median = sigma[ii] # % median of Sa and standard deviation of lnSa
            rho = baker_jayaram_correlation(iT, ip) # correlation coefficient with T_star
            
            # if (ii==13):
            #     print(rho)
            #     print(ip)
            
            icms = exp(log(sa_median) + epsilon* rho* ss_median)
            cms.append(icms) # % CMS, conditioned on T_star
            cms_out[ii,0] = ip
            cms_out[ii,1] = target_sa
            
            cms_out[ii, 2*icase+2] = icms
            # print("icase=", icase," ii=",ii, " icms=", icms)
            cms_sigma.append(ss_median* sqrt(1-rho**2)) # Conditional standard deviation of lnSa
            cms_out[ii, 2*icase+3] = ss_median* sqrt(1-rho**2)
    
        # cms= np.array(cms).reshape(1,len(inPeriodList))
        # print("cms= ", cms)
        case += "Mw"+str(iM) + "_R"+ str(iR) +"_T"+ str(iT) +"," +"sigma_Mw"+str(iM) + "_R"+ str(iR) +"_T"+ str(iT) +","
        # cms_out.append([cms])
        # cms_out = np.hstack((cms_out,cms))#, axis=1)
        # print("caseout ",cms_out)
    
    fname = outcsvfilename #'CMS_NGA_east.csv'
    with open(fname, 'w') as f:    
        # np.savetxt(f, case, delimiter=',')
        f.write(case[:-1] + '\n')
    with open(fname, 'ba') as f:
        np.savetxt(f, cms_out, delimiter=',')
    # return()
        
####################################################### 
################ USER INPUT ########################### 
####################################################### 
    
#### 
"""
compute CMS
## hazard curves for inidvidual GMM are NOT available
"""
###
# inMagnitude = [7];			# Range: 4-8.2
# inRrup = [50]; 			# range: 0-1500km
inSigmaType = "Ergodic";	# options: "Ergodic" "SingleStation"
inNsigma = 1.0; 

vs30 = 2000.0   #% site conditions, as measured by Vs30
# vs30 = 3000
# inPeriodList = np.array([0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.075,
#                 0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,
#                 1,1.5,2,3,4,5,7.5,10])#,"PGA","PGV"); # range:0.01-10sec

inPeriodList = np.array([0.01,0.02,0.03,0.05,0.075,
                0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,
                1,1.5,2,3,4,5,7.5,10])#,"PGA","PGV"); # range:0.01-10sec

# ##% target event
# # T_star = [0.1, 0.5, 1, 2, 5]; #% period at which PSHA was performed (the rest of the CMS is conditioned on this)

# ## 10,000year UHS for Lynn
# UHS_periods = [0.01,0.02,0.03,0.05,0.07,0.10,0.15,0.20,0.25,0.30,0.40,
#                0.50,0.75,1.00,1.50,2.00,3.00,4.00,5.00,7.50,10.00]
# UHS_sa = [0.4298,0.7187,0.8549,0.9594,0.8767,0.7060,0.6197,0.4846,0.3964,
#           0.3367,0.2532,0.2030,0.1367,0.1005,0.0621,0.0445,0.0256,0.0172,
#           0.0125,0.0072,0.0047]

# ## 10,000year UHS for Salem (2,000m/s)
# UHS_periods = [0.00,0.01,0.02,0.03,0.05,0.07,0.10,0.15,0.20,0.25,0.30,0.40,
#                0.50,0.75,1.00,1.50,2.00,3.00,4.00,5.00,7.50,10.00]
# UHS_sa = [0.4662,0.4662,0.7815,0.9289,1.0329,0.9452,0.7562,0.6622,
#           0.5151,0.4209,0.3558,0.2665,0.2124,0.1421,0.1039,0.0639,
#           0.0456,0.0262,0.0176,0.0127,0.0073,0.0047]

### create a list of cases
# list_mag = [5,7,5,7]
# list_R = [10,10,50,50]

# # list_mag = [7, 7, 7, 7, 7]
# # list_R = [50, 50, 50, 50, 50]

# list_Tstar = [2,2,2,2]

### for multiple UHS
uhs_folder = r"C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\Bedrock PSHA results\Salem\UHRS/"
list_rp = [10, 50, 100, 250, 475, 975, 2475, 5000, 10000]
# list_rp = [#10, 50, 100, 250, 
#            475, 975, 2475, 5000, 10000]
### read deagg results
deagg_file = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\Bedrock PSHA results/Salem_deagg_summary.csv'
deagg = np.loadtxt(deagg_file, delimiter=",", skiprows=1)
# print(deagg)

# list_Tstar0 = [0.2,1,5]

for irp in list_rp:
    list_Tstar0 = [0.2,1,5]
    uhs_file = uhs_folder+'UHS_'+str(irp)+'yr.txt'
    uhs_data = np.loadtxt(uhs_file)
    # print(uhs_data)
    # print(type(uhs_data))
    UHS_periods = uhs_data[:,0].tolist()
    UHS_sa = uhs_data[:,1].tolist()

    mask = deagg[:,0]==irp
    deagg_sub = deagg[mask,:]
    
    list_mag = []
    list_R = []
    list_Tstar = []
    for it in list_Tstar0:
        ind = np.where(deagg_sub[:,1]==it)
        print(ind)
        print(type(ind))
        print(len(ind))
        print("ind=", deagg_sub[ind,2].size)
        # print("ind=", ind.size)
        
        if deagg_sub[ind,2].size<1:
            continue
        else:
            list_Tstar.append(it)
            list_mag.append(deagg_sub[ind,2]) #[6.3, 7.0, 7.5]
            list_R.append(deagg_sub[ind,3]) #[43, 117, 530]
        # list_Tstar = deagg_sub[:,1] #[0.5, 5, 5]
    print(list_Tstar)
    print(len(list_Tstar))
    if len(list_Tstar)>0:
        # print(type(list_mag))
        # print(type(list_R))
        # print(type(UHS_sa))
        outcsvfilename = f'CMS_Salem_{irp}yr.csv'
        
        CMS_calc(list_Tstar, inPeriodList, 
                    list_mag, list_R, 
                    UHS_periods, UHS_sa,
                    outcsvfilename)
    
        ### plot cms
        ## read the CSV and plot results  
        
        N = len(list_Tstar)
        list_colors = sns.color_palette("hls", N)
        list_ltypes = ["solid", "dashed","solid", "dashed","solid", "dashed",
                        "solid", "dashed","solid", "dashed","solid", "dashed"]
        fname = outcsvfilename #'CMS_NGA_east.csv'
        # Data for plotting
        df = pd.read_csv(fname)
        # print(df.head)
        
        t = df.iloc[:,0]
        
        fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, sharex=True,
                                            figsize=(12, 6))
        
        # ax.loglog(df)
        
        ax0.loglog(df[df.columns[0]], df[df.columns[1]], "k-", label=df.columns[1])
        ax0.loglog(df[df.columns[0]], 0.75*df[df.columns[1]], "k--", 
                    linewidth=1,
                    label='target (75%)')
        
        for ii in range(0,N):
            ax0.loglog(df[df.columns[0]], df[df.columns[ii*2+2]], 
                      color = list_colors[ii], 
                      # marker='o', 
                      linestyle=list_ltypes[ii],
                      linewidth=1, #markersize=12,
                      label=df.columns[ii*2+2])
        
            
        # ax.loglog(t, df.iloc[:,3], HSV_tuples[2], label=df.columns[3])
        
        ax0.set(xlabel='Period (s)', ylabel='Sa (g)',
                title='UHS vs CMS')
        ax0.grid()
        
        legend = ax0.legend(loc='lower left', shadow=True, fontsize=12)
        # Put a nicer background color on the legend.
        # legend.get_frame().set_facecolor('C0')
        
        
        ax1.semilogx(df[df.columns[0]], df[df.columns[1]], "k-", label=df.columns[1])
        ax1.semilogx(df[df.columns[0]], 0.75*df[df.columns[1]], "k--", 
                      linewidth=1,
                      label='target (75%)')
        
        for ii in range(0,N):
            ax1.semilogx(df[df.columns[0]], df[df.columns[ii*2+2]], 
                      color = list_colors[ii], 
                      # marker='o', 
                      linestyle=list_ltypes[ii],
                      linewidth=1, #markersize=12,
                      label=df.columns[ii*2+2])
        
        ax1.set(xlabel='Period (s)', ylabel='Sa (g)',
                title='UHS vs CMS')
        ax1.grid()
        
        
        fig.savefig(f'CMS_Salem_{irp}yr.png', dpi=300)
        plt.show()
            
        
        ## hazard curves for inidvidual GMM are available 


#### plot all CMS
fig, ax0 = plt.subplots(nrows=1, ncols=1, sharex=True,
                                    figsize=(7, 6))

N = len(list_rp)
list_colors = sns.color_palette("hls", N)

icount = -1
for irp in list_rp:  
    icount +=1
    outcsvfilename = f'CMS_Salem_{irp}yr.csv'
    # cms = np.genfromtxt(outcsvfilename,delimiter=',',dtype=None, names=True)
    datain  = pd.read_csv(outcsvfilename)
    cms = datain.values
    colnames  = datain.columns.tolist()
    
    ax0.loglog(datain[datain.columns[0]], datain[datain.columns[1]], 
               color = list_colors[icount], 
                  # marker='o', 
               linestyle='dashed',
               linewidth = 0.5, #label=df.columns[1]
               )
    for kk in range(2,len(colnames),2):
        ax0.loglog(datain[datain.columns[0]], datain[datain.columns[kk]], 
                   color = list_colors[icount], #"k--", 
               linewidth =2, #label=df.columns[1]
               )
    
ax0.set(xlabel='Period (s)', ylabel='Sa (g)',
            title='CMS')
ax0.grid()


fig.savefig(f'CMS_Salem_summary.png', dpi=300)
plt.show()
    
    







