# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 22:41:37 2021

@author: FELi

## estimate SD5-75 duration 
"""


from math import log, exp

def LG14_CENA_D75(m, r, soil):
  
##  Lee, J. and Green, R.A., 2014. An empirical significant duration relationship for stable continental regions. Bulletin of earthquake engineering, 12(1), pp.217-235.
## soil = 0 for rock sites (e.g. site class A and B) and soil=1 for soil sites  

    param = "SD5-75, H CENA"
    c1 = 0
    c2 = 2.23
    c3 = 0.10
    s1 = -0.72
    s2 = -0.19
    s3 = -0.014
    tauln = 0.46
    sigmaln = 0.35
    sigmaTln = 0.58
    # stringsAsFactors = F
  
    lnSD = log(c1 + c2*exp(m-6) + c3*r + (s1 + s2*(m-6) + s3*r)*(soil))
    
    median = exp(lnSD)
    v16th = exp(lnSD-sigmaTln)
    v84th = exp(lnSD+sigmaTln)
    
    return(m, r, median,  sigmaTln, v16th, v84th)
  
M = 5. #[5.]
R = 10. #[10.]
soil = 0
out = LG14_CENA_D75(M, R, soil)
print(out)