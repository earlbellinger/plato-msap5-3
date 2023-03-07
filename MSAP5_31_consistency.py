#!usr/bin/env python 
#### PLATO MSAP5-31 Consistency checks 
#### Author: Earl Patrick Bellinger ( bellinger@phys.au.dk ) 
#### Max Planck Institute for Astrophysics, Garching, Germany 
#### Stellar Astrophysics Centre, Aarhus, Denmark 

import numpy as np
from scipy.stats import tukey_hsd

def consistency(IDP_MASS_SEISMIC=None, 
                IDP_MASS_GRANULATION=None, 
                IDP_MASS_GRANULATION_CGBM=None, 
                IDP_MASS_RHO_TRANSIT_CGBM=None, 
                IDP_RADIUS_SEISMIC=None, 
                IDP_RADIUS_GRANULATION_CGBM=None, 
                IDP_RADIUS_RHO_TRANSIT=None, 
                IDP_RADIUS_RHO_TRANSIT_CGBM=None, 
                IDP_AGE_SEISMIC=None, 
                IDP_AGE_GYRO=None, 
                IDP_AGE_ACTIVITY=None, 
                IDP_AGE_GRANULATION_CGBM=None, 
                IDP_AGE_RHO_TRANSIT_CGBM=None):
    """
    Determines the consistency in mass, radius, and age across several measurement methods on a star-by-star basis using the Tukey Honestly Significantly Difference (HSD) test. 
    
    Parameters
    ----------
    IDP_MASS_... : array
    	Samples from the posterior mass distribution from each method supplied, in solar units 
    IDP_RADIUS_... : array
    	Samples from the posterior radius distribution from each method supplied, in solar units 
    IDP_AGE_... : array
    	Samples from the posterior age distribution from each method supplied, in Gyr 
    """
    
    Ms = [IDP_MASS_SEISMIC, 
          IDP_MASS_GRANULATION, 
          IDP_MASS_GRANULATION_CGBM, 
          IDP_MASS_RHO_TRANSIT_CGBM] 
    
    Rs = [IDP_RADIUS_SEISMIC, 
          IDP_RADIUS_GRANULATION_CGBM, 
          IDP_RADIUS_RHO_TRANSIT, 
          IDP_RADIUS_RHO_TRANSIT_CGBM]
    
    As = [IDP_AGE_SEISMIC, 
          IDP_AGE_GYRO, 
          IDP_AGE_ACTIVITY, 
          IDP_AGE_GRANULATION_CGBM, 
          IDP_AGE_RHO_TRANSIT_CGBM]
    
    # remove measurements that are None 
    Ms = [M for M in Ms if M is not None and len(M)>1]
    Rs = [R for R in Rs if R is not None and len(R)>1]
    As = [A for A in As if A is not None and len(A)>1]
    
    #print(Rs)
    
    flags = []
    for xs in [Ms, Rs, As]:
        #print(1)
        
        #print(xs)
        
        # case: no measurements (all None)
        if xs == []:
            flags += [False]
            continue 
        
        # perform test 
        test = tukey_hsd(*xs)
        #print(test.pvalue)
        flags += [np.all(test.pvalue >= 0.05)]
    
    return flags 
