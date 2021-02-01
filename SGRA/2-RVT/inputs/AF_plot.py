# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 12:01:28 2021

@author: FELi
plot a selected param, e.g. AF, over bedrock motions

"""

import numpy as np
import matplotlib.pyplot as plt
# from math import median
from statistics import median
# from math import exp
import matplotlib.backends.backend_pdf
from scipy.optimize import curve_fit 
from math import exp, log, sqrt

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

### fitting AF curves
###    ln(AF) = f1+f2*ln((sa_rock+f3)/f3)
folder_parent = r'C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\20368492 National Grid\SGRA/'
rvt_output_folder = folder_parent+'2-RVT/'

cms_initial = 'Salem_cms_initial.csv'
param = "Amplification Factor"
param_file = 'RVT1-0.00 ft (Outcrop (2A)) from Bedrock (Outcrop (2A))-specRatio.csv'
coeff = 1.0
# param = "CSR"
# param_file = 'profile-stressRatio.csv'
# coeff = 0.65

param_data = np.genfromtxt(rvt_output_folder+param_file, delimiter=',', dtype=None,#names=True, 
                        skip_header=3)
cmsdata  = np.genfromtxt(rvt_output_folder+'inputs/'+cms_initial, delimiter=',', dtype=None, names=True)

# print(len(cmsdata))
no_periods = len(cmsdata)
no_cms = len(cmsdata[0].tolist())-1

# print(type(param_data))
# print(type(param_data[0]))
# print(param_data[0])

param_data = np.array(param_data.tolist())

no_params = len(param_data[0,:])-3

no_sims = no_params/no_cms
# print(AF_data[:,0])
pdf = matplotlib.backends.backend_pdf.PdfPages(rvt_output_folder+f"{param}_output.pdf")

ii = 0
curfit_params = []
AF_median_cms = []
for ii in range(no_periods):

    temp = cmsdata[ii].tolist()
    ip = temp[0]
    print(ip)
    # print(type(temp.tolist()))
    sa_bedrock = np.asarray(temp[1:])

    param_periods = param_data[:,0]
    # print(AF_periods)
    
    idx,param_p = find_nearest(param_periods, ip)
    
    fig, ax0 = plt.subplots(nrows=1, ncols=1, sharex=True,
                                                figsize=(7, 6))
    median_out = [ip]
    sa_bedrock_all = []
    param_all = []
    
    for jj in range(no_cms):
        temp_idx = range(jj+1,no_params+1, no_cms)
        param_datasub = param_data[idx,temp_idx]*coeff
        
        median_out.append(median(param_datasub))
        
        temp_x = np.repeat(sa_bedrock[jj],no_sims)

        ax0.semilogx(temp_x, param_datasub, "ko")
        sa_bedrock_all.append(temp_x)
        param_all.append(param_datasub)
    
    AF_median_cms.append(median_out)
    
    sa_bedrock_all = np.asarray(sa_bedrock_all).flatten()
    param_all = np.asarray(param_all).flatten()
    # print(type(sa_bedrock_all))
    # print(type(param_all))
    # print(param_all)
    p = sa_bedrock_all.argsort()
    
    ## polynomial fit
    x=  sa_bedrock_all[p]
    y = np.asarray(param_all)[p]
    # ax0.semilogx(x,y, "k-")

    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    xp = np.logspace(-4, 0, 50)
    ax0.semilogx(xp,p(xp), "r-", label='polynomial-fit')
    
    ## curve fit
    # x= x.tolist()
    
    def obj(x, f1, f2, f3):
        
        return np.exp(f1+f2*np.log((x+f3)/f3))
    
    param_value, param_cov = curve_fit(obj, x, y) 
    # print(type(param_value))
    param_value = np.round(param_value,3)
    param_cov = np.round(param_cov,3)
    
    # print("funcion coefficients:", param_value) 
    # print("Covariance of coefficients:", param_cov) 
    ax0.text(0.0001, 3.5, "f1="+ str( param_value[0])+
             "; f2="+ str( param_value[1])+"; f3="+ str( param_value[2]),
             fontsize=10)
    # ax0.text(0.0001, 3.0, "COV:"+ str(param_cov),
    #          fontsize=10)
    # ans = (param[0]*(np.sin(param[1]*x))) 
    ans = np.exp(param_value[0]+param_value[1]*np.log((x+param_value[2])/param_value[2]))
    # ax0.plot(x, y, 'o', color ='red', label ="data") 
    SD_res = round(sqrt(sum((y-ans)**2)/(len(x)-3)),3)
    ax0.text(0.0001, 3.0, "SD: "+ str(SD_res),
             fontsize=10)
    
    curfit_params.append([ip, param_value[0], param_value[1],
                           param_value[2], SD_res, sqrt(param_cov[0,0]),
                          sqrt(param_cov[1,1]),
                          sqrt(param_cov[2,2])])
    ax0.plot(x, ans, '--', color ='blue', 
             linewidth=2,
             label ="Stewart-fit") 
    ax0.legend() 

    ####----
    ax0.set_ylim([0,4])
    ax0.set(xlabel='SA_bedrock (g)', ylabel=f'{param}',
                    title=f'{ip}')
    ax0.grid()
    
    # fig.savefig(rvt_output_folder+f'{param}_Salem_{ii}_{ip}s.png', dpi=300)
    plt.show()
    
    pdf.savefig( fig )
    
pdf.close()     

fname = rvt_output_folder+'fit_aprams.csv' #'CMS_NGA_east.csv'
with open(fname, 'w') as f:    
    # np.savetxt(f, case, delimiter=',')
    f.write('period_s,f1,f2,f3,sd_res,cov1,cov2,cov3' + '\n')
with open(fname, 'ba') as f:
    np.savetxt(f, curfit_params, delimiter=',')

fname = rvt_output_folder+'AF_median_cms.csv' #'CMS_NGA_east.csv'
temp_str = 'period_s'
for ii in range(no_cms):
    temp_str += ','+str(ii+1)

with open(fname, 'w') as f:    
    # np.savetxt(f, case, delimiter=',')
    f.write(temp_str + '\n')
with open(fname, 'ba') as f:
    np.savetxt(f, AF_median_cms, delimiter=',')

