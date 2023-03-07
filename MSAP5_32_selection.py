#!usr/bin/env python 
#### PLATO MSAP5-32 Selection
#### Author: Earl Patrick Bellinger ( bellinger@phys.au.dk ) 
#### Max Planck Institute for Astrophysics, Garching, Germany 
#### Stellar Astrophysics Centre, Aarhus, Denmark 

import numpy as np
from MSAP5_31_consistency import consistency

M_names = ['IDP_MASS_SEISMIC', 
           'IDP_MASS_GRANULATION', 
           'IDP_MASS_GRANULATION_CGBM', 
           'IDP_MASS_RHO_TRANSIT_CGBM'] 

R_names = ['IDP_RADIUS_SEISMIC', 
           'IDP_RADIUS_GRANULATION_CGBM', 
           'IDP_RADIUS_RHO_TRANSIT', 
           'IDP_RADIUS_RHO_TRANSIT_CGBM']

A_names = ['IDP_AGE_SEISMIC', 
           'IDP_AGE_GYRO', 
           'IDP_AGE_ACTIVITY', 
           'IDP_AGE_GRANULATION_CGBM', 
           'IDP_AGE_RHO_TRANSIT_CGBM']

def selection(IDP_MASS_SEISMIC=None, 
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
    Selects from a priority list the mass, radius, and age to present as the final measurement. 
    
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
    
    consistent = consistency(*Ms, *Rs, *As)
    
    M_result = None
    if consistent[0]:
        for M_name, M in zip(M_names, Ms):
            if M is None or len(M)<=1:
                continue
            M_result = M_name
            break 
    
    R_result = None
    if consistent[1]:
        for R_name, R in zip(R_names, Rs):
            if R is None or len(R)<=1:
                continue
            R_result = R_name
            break 
    
    A_result = None
    if consistent[2]:
        for A_name, A in zip(A_names, As):
            if A is None or len(A)<=1:
                continue
            A_result = A_name
            break 
    
    return [M_result, R_result, A_result]
