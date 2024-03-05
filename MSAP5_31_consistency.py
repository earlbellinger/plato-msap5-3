#!usr/bin/env python 
#### PLATO MSAP5-31 Consistency checks 
#### Author: Earl Patrick Bellinger ( earl.bellinger@yale.edu ) 
#### Department of Astronomy, Yale University

import numpy as np
from scipy.stats import tukey_hsd

def consistency(IDP_SAS_MASS_GRID_MIXED=None, 
                IDP_SAS_MASS_GRID_SURF_IND=None, 
                IDP_SAS_MASS_FREQS=None, 
                IDP_SAS_MASS_SCALING_GRIDS=None, 
                IDP_SAS_MASS_SCALING_ONLY=None, 
                IDP_SAS_MASS_GRANULATION=None, 
                IDP_SAS_MASS_GRANULATION_CGBM=None, 
                IDP_SAS_RADIUS_GRID_MIXED=None, 
                IDP_SAS_RADIUS_GRID_SURF_IND=None, 
                IDP_SAS_RADIUS_FREQS=None, 
                IDP_SAS_RADIUS_SCALING_GRIDS=None, 
                IDP_SAS_RADIUS_SCALING_ONLY=None, 
                IDP_SAS_RADIUS_GRANULATION_CGBM=None, 
                IDP_SAS_AGE_GRID_MIXED=None, 
                IDP_SAS_AGE_GRID_SURF_IND=None, 
                IDP_SAS_AGE_FREQS=None, 
                IDP_SAS_AGE_SCALING_GRIDS=None, 
                IDP_SAS_AGE_GYRO=None, 
                IDP_SAS_AGE_ACTIVITY=None, 
                IDP_SAS_AGE_GRANULATION_CGBM=None):
    """
    Determines the consistency in mass, radius, and age across several measurement methods on a star-by-star basis using the Tukey Honestly Significantly Difference (HSD) test. 
    
    Parameters
    ----------
    IDP_SAS_... : array
        Samples from the posterior distribution from each method supplied, in respective units (mass in solar units, radius in solar units, age in Gyr)
    """

    Ms = [IDP_SAS_MASS_GRID_MIXED, 
          IDP_SAS_MASS_GRID_SURF_IND, 
          IDP_SAS_MASS_FREQS, 
          IDP_SAS_MASS_SCALING_GRIDS, 
          IDP_SAS_MASS_SCALING_ONLY, 
          IDP_SAS_MASS_GRANULATION, 
          IDP_SAS_MASS_GRANULATION_CGBM]

    Rs = [IDP_SAS_RADIUS_GRID_MIXED, 
          IDP_SAS_RADIUS_GRID_SURF_IND, 
          IDP_SAS_RADIUS_FREQS, 
          IDP_SAS_RADIUS_SCALING_GRIDS, 
          IDP_SAS_RADIUS_SCALING_ONLY, 
          IDP_SAS_RADIUS_GRANULATION_CGBM]

    As = [IDP_SAS_AGE_GRID_MIXED, 
          IDP_SAS_AGE_GRID_SURF_IND, 
          IDP_SAS_AGE_FREQS, 
          IDP_SAS_AGE_SCALING_GRIDS, 
          IDP_SAS_AGE_GYRO, 
          IDP_SAS_AGE_ACTIVITY, 
          IDP_SAS_AGE_GRANULATION_CGBM]
    
    # remove measurements that are None 
    Ms = [M for M in Ms if M is not None and len(M)>1]
    Rs = [R for R in Rs if R is not None and len(R)>1]
    As = [A for A in As if A is not None and len(A)>1]
    
    flags = []
    for xs in [Ms, Rs, As]:
        # case: no measurements (all None)
        if xs == []:
            flags += [False]
            continue 
        
        # perform test 
        test = tukey_hsd(*xs)
        flags += [np.all(test.pvalue >= 0.05)]
    
    return flags 
