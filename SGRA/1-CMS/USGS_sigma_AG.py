# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 23:30:49 2021

@author: FELi
"""
import numpy as np
from scipy.interpolate import interp1d
import math



def sigmaEpri(Mw, ip):
    
    ## # T,     tau_M5, phi_M5, tau_M6, phi_M6, tau_M7, phi_M7
    coeff_sigma_usgs2018_epri = np.array([  
      [0.01,  0.4320, 0.6269, 0.3779, 0.5168, 0.3525, 0.5039],
      [0.02,  0.4710, 0.6682, 0.4385, 0.5588, 0.4138, 0.5462],
      [0.03,  0.4710, 0.6682, 0.4385, 0.5588, 0.4138, 0.5462],
    # 0.04,  0.4710, 0.6682, 0.4385, 0.5588, 0.4138, 0.5462
      [0.05,  0.4710, 0.6682, 0.4385, 0.5588, 0.4138, 0.5462],
      [0.075, 0.4710, 0.6682, 0.4385, 0.5588, 0.4138, 0.5462],
      [0.1,   0.4710, 0.6682, 0.4385, 0.5588, 0.4138, 0.5462],
      [0.15,  0.4433, 0.6693, 0.4130, 0.5631, 0.3886, 0.5506],
      [0.2,   0.4216, 0.6691, 0.3822, 0.5689, 0.3579, 0.5566],
      [0.25,  0.4150, 0.6646, 0.3669, 0.5717, 0.3427, 0.5597],
      [0.3,   0.4106, 0.6623, 0.3543, 0.5846, 0.3302, 0.5727],
      [0.4,   0.4088, 0.6562, 0.3416, 0.5997, 0.3176, 0.5882],
      [0.5,   0.4175, 0.6526, 0.3456, 0.6125, 0.3217, 0.6015],
      [0.75,  0.4439, 0.6375, 0.3732, 0.6271, 0.3494, 0.6187],
      [1.0,   0.4620, 0.6219, 0.3887, 0.6283, 0.3650, 0.6227],
      [1.5,   0.4774, 0.5957, 0.4055, 0.6198, 0.3819, 0.6187],
      [2.0,   0.4809, 0.5860, 0.4098, 0.6167, 0.3863, 0.6167],
      [3.0,   0.4862, 0.5813, 0.4186, 0.6098, 0.3952, 0.6098],
      [4.0,   0.4904, 0.5726, 0.4144, 0.6003, 0.3910, 0.6003],
      [5.0,   0.4899, 0.5651, 0.4182, 0.5986, 0.3949, 0.5986],
      [7.5,   0.4803, 0.5502, 0.4067, 0.5982, 0.3835, 0.5982],
      [10.0,  0.4666, 0.5389, 0.3993, 0.5885, 0.3761, 0.5885],
     # [0.00,   0.4320, 0.6269, 0.3779, 0.5168, 0.3525, 0.5039],
    # PGV,   0.3925, 0.5979, 0.3612, 0.5218, 0.3502, 0.5090]
    ])


   ## tau_M5, phi_M5, tau_M6, phi_M6, tau_M7, phi_M7
    T = coeff_sigma_usgs2018_epri[:,0]
    tau_M5_all = coeff_sigma_usgs2018_epri[:,1]
    phi_M5_all = coeff_sigma_usgs2018_epri[:,2]
    tau_M6_all = coeff_sigma_usgs2018_epri[:,3]
    phi_M6_all = coeff_sigma_usgs2018_epri[:,4]
    tau_M7_all = coeff_sigma_usgs2018_epri[:,5]
    phi_M7_all = coeff_sigma_usgs2018_epri[:,6]
    
    if (Mw <= 5.0):
        tau_all = tau_M5_all
        phi_all = phi_M5_all
    elif (Mw <= 7.0):
        tau_all = interp1d(np.array([5.0, 6.0, 7.0]),np.transpose(np.array([tau_M5_all, tau_M6_all, tau_M7_all])))(Mw)
        phi_all = interp1d(np.array([5.0, 6.0, 7.0]),np.transpose(np.array([phi_M5_all, phi_M6_all, phi_M7_all])))(Mw)
    else:
        tau_all = tau_M7_all
        phi_all = phi_M7_all
    # print("ip=",ip, ", tau=", tau, ",phi=",phi)
    tau = interp1d(T, tau_all)(ip)
    phi = interp1d(T, phi_all)(ip)
    return np.hypot(phi, tau)

def sigmaPanel(Mw, vs30, ip):
    
    
    ### #      T,    t1,  t2,  t3,  t4, ss_a, ss_b,  s2s1,  s2s2
    coeff_sigma_usgs2018_panel = np.array([
     
      [0.01,  0.4436, 0.4169, 0.3736, 0.3415, 0.5423, 0.3439, 0.533, 0.566],
      [0.02,  0.4436, 0.4169, 0.3736, 0.3415, 0.5410, 0.3438, 0.537, 0.577],
      [0.03,  0.4436, 0.4169, 0.3736, 0.3415, 0.5397, 0.3437, 0.542, 0.598],
    # [0.04,  0.4436, 0.4169, 0.3736, 0.3415, 0.5382, 0.3436, 0.562, 0.638],
      [0.05,  0.4436, 0.4169, 0.3736, 0.3415, 0.5371, 0.3435, 0.583, 0.653],
      [0.075, 0.4436, 0.4169, 0.3736, 0.3415, 0.5339, 0.3433, 0.619, 0.633],
      [0.1,   0.4436, 0.4169, 0.3736, 0.3415, 0.5308, 0.3431, 0.623, 0.590],
      [0.15,  0.4436, 0.4169, 0.3736, 0.3415, 0.5247, 0.3466, 0.603, 0.532],
      [0.2,   0.4436, 0.4169, 0.3736, 0.3415, 0.5189, 0.3585, 0.578, 0.461],
      [0.25,  0.4436, 0.4169, 0.3736, 0.3415, 0.5132, 0.3694, 0.554, 0.396],
      [0.3,   0.4436, 0.4169, 0.3736, 0.3415, 0.5077, 0.3808, 0.527, 0.373],
      [0.4,   0.4436, 0.4169, 0.3736, 0.3415, 0.4973, 0.4004, 0.491, 0.339],
      [0.5,   0.4436, 0.4169, 0.3736, 0.3415, 0.4875, 0.4109, 0.472, 0.305],
      [0.75,  0.4436, 0.4169, 0.3736, 0.3415, 0.4658, 0.4218, 0.432, 0.273],
      [1.0,   0.4436, 0.4169, 0.3736, 0.3415, 0.4475, 0.4201, 0.431, 0.257],
      [1.5,   0.4436, 0.4169, 0.3736, 0.3415, 0.4188, 0.4097, 0.424, 0.247],
      [2.0,   0.4436, 0.4169, 0.3736, 0.3415, 0.3984, 0.3986, 0.423, 0.239],
      [3.0,   0.4436, 0.4169, 0.3736, 0.3415, 0.3733, 0.3734, 0.418, 0.230],
      [4.0,   0.4436, 0.4169, 0.3736, 0.3415, 0.3604, 0.3604, 0.412, 0.221],
      [5.0,   0.4436, 0.4169, 0.3736, 0.3415, 0.3538, 0.3537, 0.404, 0.214],
      [7.5,   0.4436, 0.4169, 0.3736, 0.3415, 0.3482, 0.3481, 0.378, 0.201],
      [10.0,  0.4436, 0.4169, 0.3736, 0.3415, 0.3472, 0.3471, 0.319, 0.193],
    #  [0.00,   0.4436, 0.4169, 0.3736, 0.3415, 0.5423, 0.3439, 0.533, 0.566],
    # [PGV,   0.3633, 0.3532, 0.3340, 0.3136, 0.4985, 0.3548, 0.487, 0.458]
      ])

    
    #### T,    t1,  t2,  t3,  t4, ss_a, ss_b,  s2s1,  s2s2
    T = coeff_sigma_usgs2018_panel[:,0]
    tau1_all = coeff_sigma_usgs2018_panel[:,1]
    tau2_all = coeff_sigma_usgs2018_panel[:,2]
    tau3_all = coeff_sigma_usgs2018_panel[:,3]
    tau4_all = coeff_sigma_usgs2018_panel[:,4]
    ss_a_all = coeff_sigma_usgs2018_panel[:,5]
    ss_b_all = coeff_sigma_usgs2018_panel[:,6]
    phi_s2s1_all = coeff_sigma_usgs2018_panel[:,7]
    phi_s2s2_all = coeff_sigma_usgs2018_panel[:,8]
        
    tau_all = tau_calc(Mw, tau1_all, tau2_all, tau3_all, tau4_all)
    phi_ss_all = phi_ss_calc(Mw, ss_a_all, ss_b_all);

    ### φ_s2s model; single branch; Stewart et al. (2019) */
    phi_s2s_all = phi_s2s_calc(vs30, phi_s2s1_all, phi_s2s2_all);
    # print("tau=", tau, ", phi_ss=" , phi_ss, ",phi_s2s=", phi_s2s)

    tau = interp1d(T, tau_all)(ip)
    phi_ss = interp1d(T, phi_ss_all)(ip)
    phi_s2s = interp1d(T, phi_s2s_all)(ip)

    return np.sqrt(tau**2+phi_ss**2+ phi_s2s**2);

def tau_calc(Mw, tau1, tau2, tau3, tau4):
    if (Mw <= 4.5):
      return tau1;
    elif (Mw <= 5.0):
      return tau1 + (tau2 - tau1) * (Mw - 4.5) / 0.5;
    elif (Mw <= 5.5):
      return tau2 + (tau3 - tau2) * (Mw - 5.0) / 0.5;
    elif (Mw <= 6.5):
      return tau3 + (tau4 - tau3) * (Mw - 5.5);
    
    return tau4;


def phi_ss_calc(Mw, ss_a, ss_b):
    if (Mw <= 5.0):
      return ss_a;
    elif (Mw <= 6.5):
      return ss_a + (Mw - 5.0) * (ss_b - ss_a) / 1.5;
    
    return ss_b;

def phi_s2s_calc(vs30, phi_s2s1, phi_s2s2):
    V_phi1 = 1200.0;
    V_phi2 = 1500.0;
    if (vs30 < V_phi1):
      return phi_s2s1;
    elif (vs30 < V_phi2):
      return phi_s2s1 - ((phi_s2s1 - phi_s2s2) / (V_phi2 - V_phi1)) * (vs30 - V_phi1);
    
    return phi_s2s2;


def sigmatotal(Mw, vs30,
               period_list):
    SIGMA_LTC_WTS = [ 0.2, 0.8 ]; ##for panel and epri
    s1 = sigmaPanel(Mw, vs30, period_list)
    s2 = sigmaEpri(Mw, period_list)
    sigma = s1*SIGMA_LTC_WTS[0] + s2*SIGMA_LTC_WTS[1]
        
    # print("ip=",ip, "s1=",s1, ", s2=",s2)
    return sigma


if __name__ == "__main__":
    Mw = 5
    vs30 = 2000
    period_list = [0.01, 0.1,1,5]
    period_list = (0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.075,
                    0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,
                    1,1.5,2,3,4,5,7.5,10)
    sigma = sigmatotal( 
                    Mw, vs30,
                    period_list)     
    print("sigma= ", sigma)



