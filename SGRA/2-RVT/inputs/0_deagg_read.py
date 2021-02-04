# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 11:04:45 2021

@author: FELi
"""
import numpy as np

# period = 
# deagg_folder = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\Bedrock PSHA results\Salem\deagg/'
# deagg_file = deagg_folder+'10yr/summary.xlsx'


# x= np.array([10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,50,50,50,50,50,50,50,50,50,50,50,50,50,50,70,70,70,70,70,70,70,70,70,70,90,90,90,90,90,90,90,90,90,110,110,110,110,110,110,110,110])
# y = np.array([4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9])
# x = (x-10)/20
# x = x.astype(int)


# y= (y-4.7)/0.2
# y = np.rint(y)
# # print(round(y,2))
# y = y.astype(np.int)

# print(y)

# z = [9.56,10.04,10.6,11.08,11.56,11.94,12.3,11.57,12.47,12.42,12.22,11.93,12.69,13.22,13.1,13.15,13.18,23.16,23.94,24.79,25.54,26.15,26.51,26.85,26.68,29.05,28.8,28.38,27.76,29.45,30.05,29.36,29.84,30.74,42.79,43.02,43.55,44.73,48.61,48.8,48.39,47.67,47.48,49.48,49.17,49.36,50.04,50.04,64.83,67.68,68.34,67.71,69.12,70.03,69.13,68.31,69.69,68.06,83.29,86.65,86.89,88.57,89.59,88.35,88.53,88.45,89.05,102.64,106.06,108.52,109.54,108.27,109.28,108.55,109.46]
# # xyz = np.array((x,y,z)).T
# # print(xyz[:,1])
# # # arr = np.zeros((6,11))
# # # yx = zip(y,x)

# z_array=np.zeros((17,6),float); 
# # z_array=np.zeros((6,17),float); 
# z_array[y,x]=z


# print (z_array)


import glob
import pandas as pd
import numpy as np
from math import log, exp
import matplotlib.pyplot as plt

# from EstimateDuration import LG14_CENA_D75

folder = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\Bedrock PSHA results\Salem\deagg/'
list_rp = [10, 50, 100,250,475,975,2475,5000,10000]
# list_folders = glob.glob(folder+'*/')
# print(list_folders)

# ifolder = list_folders[0]
# deagg_out = []


M4matrix = np.linspace(4.7,7.9,17).round(1)
R4matrix = np.linspace(10,990,50)

ip = 0.01 #'PGA'

## read hazard curves
hc_folder = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\Bedrock PSHA results\Salem\hazard curves/'

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
print('aep_bedrock=',aep_bedrock)
out_sa = []
out_con =[]
for irp in  list_rp:
    ifolder =folder + f'{irp}yr/' # list_folders:
    print(ifolder)
       
    # rp = ifolder.split('\\')[-2]
    yr = float(irp)
    print(log(1/yr))
    print('1/yr=',1/yr)
    temp = np.interp(log(1/yr), np.log(aep_bedrock)[::-1], np.log(sa_bedrock)[::-1])
    sa = exp(temp)
    print('temp=', temp,"  sa=", sa)
    
    out_sa.append(sa)
    
    excelfile = ifolder+'summary.xlsx'
    xl = pd.ExcelFile(excelfile)
    # print(xl.sheet_names)
    
    sht_names = xl.sheet_names
    if ip ==0.01:
        isht = 'PGA'
    else:
        isht = f'{ip}s'
    if isht in sht_names:
        
    # for isht in sht_names:
    #     if isht=='PGA':
    #         period=0.01
    #     else:
    #         period = float(isht.replace('s',''))
        df = xl.parse(isht)
        if len(df.index)<1:
            continue
        else:        
            print("isht=", isht)
            # print(df.iloc[5:20,0:6])
            temp = df.iloc[5,0]
            mw = float(temp.split(':')[1])
            temp = df.iloc[6,0]
            r = float(temp.split(':')[1].split(" ")[2])
            
            # print("mw=",mw, ", r=", r, " ,", type(mw))   
            
            list_m =  df.iloc[17:,2]
            list_r =  df.iloc[17:,0]
            list_con =  df.iloc[17:,5]
            
            cons = list_con.values
            con_out = []
            for ii in cons:
                # print(ii, ', ', type(ii))
                
                if isinstance(ii, str):
                    con_out.append(0.0)
                else:
                    con_out.append(ii)
            cons = con_out
            ms = list_m.values.astype(np.float)
            rs = list_r.values.astype(np.float)
            
            ## get the matrix 
            # print('rs=',rs)
            y = np.rint((rs-10)/20)
            y = y.astype(int)
            
            # print('ms=', ms)
            x= np.rint((ms-4.7)/0.2)
            x = x.astype(np.int)
            
            # print(y)
            z_array=np.zeros((50,17),float); 
            # z_array=np.zeros((6,17),float); 
            z_array[y,x]=cons
            z_array[z_array == 0] = 'nan'
            
            con1D = z_array.flatten()
            
#     # yr
    # temp = (con1D*1./yr)
    # out_aep.append(temp.tolist())
    out_con.append(con1D.tolist())


out_con = np.array(out_con)
out_aep = []
print('out_sa=', out_sa)
print('out_aep[:,1]=',out_con[:,1])

aeps = 1./np.array(list_rp)
print('aeps=',aeps)
import seaborn as sns
N = 17
temp = sns.color_palette("hls", N)
list_colors = temp

list_ltypes = ["solid","dashed","dotted",
               "solid","dashed","dotted",
               "solid","dashed","dotted",
               "solid","dashed","dotted",
               "solid","dashed","dotted",
               "solid","dashed","dotted",] 

fit_folder = r'C:\Users\FEli\OneDrive - Golder Associates\Projects\2020\20368492 National Grid\SGRA\2-RVT/'

icount = -1
ir = R4matrix[1]
for ir in R4matrix:
    fig, ax0 = plt.subplots(nrows=1, ncols=1, sharex=True,
                                                figsize=(7, 6))
    ii=-1
    for im in M4matrix:    
        icount += 1
        ii+=1
        
        con = np.interp(aep_bedrock, aeps[::-1], out_con[:,icount][::-1])
        aep = (con*0.01*aep_bedrock)
        out_aep.append([float(ir),float(im)]+aep.tolist())
        
        if ii==0:
            print('aeps[::-1]=', aeps[::-1])
            
            print('out_con[:,icount][::-1]=', out_con[:,icount][::-1])
            
            print('con=', con)
            print('aep_bedrock=', aep_bedrock)
            
            
            print('out_aep=',aep)
        ax0.loglog(sa_bedrock, aep,
                         color =  list_colors[ii], 
                         # marker='o', 
                         linestyle=list_ltypes[ii],
                         linewidth=1, #markersize=12,
                         label= f'Mw{im}')
        
    ax0.legend() 
    
    ####----
    # ax0.set_ylim([0,1.2])
    ax0.set(xlabel='sa (g)', ylabel=f'AEP',
                    title=f'Hazard curves ({ir}km)')
    ax0.grid()
    
    fig.savefig(fit_folder+f'Salem_hc_deagg_{ir}km.png', dpi=300)
    plt.show() 


print(out_aep[1])

# deagg_out = np.array(deagg_out)
header = 'R_km,Mw'
for isa in sa_bedrock:
    isa = round(isa,5)
    header += f',sa{isa}'
fname = fit_folder+f'Salem_hc_deagg_{ip}.csv' #'CMS_NGA_east.csv'
with open(fname, 'w') as f:    
    f.write(header+'\n')
with open(fname, 'ba') as f:
    np.savetxt(f, out_aep, delimiter=',')

