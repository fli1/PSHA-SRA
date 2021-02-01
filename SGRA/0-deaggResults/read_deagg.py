# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:53:29 2021

@author: FELi
read deagg results
"""

import glob
import pandas as pd
import numpy as np
from EstimateDuration import LG14_CENA_D75

folder = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\Bedrock PSHA results\Salem\deagg/'
list_folders = glob.glob(folder+'*/')
# print(list_folders)

ifolder = list_folders[0]
deagg_out = []

## for Lee and Green D75 estimate, require this parameter:
soil_index = 0  # 0 = rock; 1 = soil

for ifolder in list_folders:
    print(ifolder)
       
    rp = ifolder.split('\\')[-2]
    yr = float(rp.replace("yr",""))
    excelfile = ifolder+'summary.xlsx'
    xl = pd.ExcelFile(excelfile)
    # print(xl.sheet_names)
    
    sht_names = xl.sheet_names
    
    for isht in sht_names:
        if isht=='PGA':
            period=0.01
        else:
            period = float(isht.replace('s',''))
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
            
            print("mw=",mw, ", r=", r, " ,", type(mw))   
            
            list_m =  df.iloc[17:,2]
            list_r =  df.iloc[17:,0]
            list_con =  df.iloc[17:,5]
            
            cons = list_con.values
            ms = list_m.values
            rs = list_r.values
            con_out = []
            for ii in cons:
                # print(ii, ', ', type(ii))
                
                if isinstance(ii, str):
                    con_out.append(0.0)
                else:
                    con_out.append(ii)
            
            con_out = np.where(rs>400.0, 0.0, con_out)
            # print(con_out)
            # print('len=', len(con_out))
            mw_mean = sum(con_out*ms)/sum(con_out)
            r_mean = sum(con_out*rs)/sum(con_out)
            print(mw_mean, " ,", r_mean, ",", type(mw_mean))
            if(r_mean/r<0.85):
                _, _, d75, d75_sigma,_, _ = LG14_CENA_D75(mw_mean, r_mean, soil_index)
                out = [yr, period, mw_mean, r_mean, d75]
            else:
                _, _, d75, d75_sigma,_, _ = LG14_CENA_D75(mw, r, soil_index)
                out = [yr, period, mw, r, d75]
            deagg_out.append(out)

# deagg_out = np.array(deagg_out)
fname = 'Salem_deagg_summary.csv' #'CMS_NGA_east.csv'
with open(fname, 'w') as f:    
    f.write('return_period_yr,period_s,mw,distance_km, D75_s\n')
with open(fname, 'ba') as f:
    np.savetxt(f, deagg_out, delimiter=',')

