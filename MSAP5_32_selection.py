#!usr/bin/env python 
#### PLATO MSAP5-32 Selection
#### Author: Earl Patrick Bellinger ( earl.bellinger@yale.edu ) 
#### Department of Astronomy, Yale University

import numpy as np
from MSAP5_31_consistency import consistency
from MSAP5 import M_names, R_names, A_names

def selection(IDP_SAS_MASS_GRID_MIXED=None, 
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
              IDP_SAS_AGE_GRANULATION_CGBM=None,
              IDP_SAS_CONSISTENCY_FLAGS=[False, False, False]):
    """
    Selects from a priority list the mass, radius, and age to present as the final measurement. 
    
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

    consistent = consistency(IDP_SAS_MASS_GRID_MIXED, 
                             IDP_SAS_MASS_GRID_SURF_IND, 
                             IDP_SAS_MASS_FREQS, 
                             IDP_SAS_MASS_SCALING_GRIDS, 
                             IDP_SAS_MASS_SCALING_ONLY, 
                             IDP_SAS_MASS_GRANULATION, 
                             IDP_SAS_MASS_GRANULATION_CGBM, 
                             IDP_SAS_RADIUS_GRID_MIXED, 
                             IDP_SAS_RADIUS_GRID_SURF_IND, 
                             IDP_SAS_RADIUS_FREQS, 
                             IDP_SAS_RADIUS_SCALING_GRIDS, 
                             IDP_SAS_RADIUS_SCALING_ONLY, 
                             IDP_SAS_RADIUS_GRANULATION_CGBM, 
                             IDP_SAS_AGE_GRID_MIXED, 
                             IDP_SAS_AGE_GRID_SURF_IND, 
                             IDP_SAS_AGE_FREQS, 
                             IDP_SAS_AGE_SCALING_GRIDS, 
                             IDP_SAS_AGE_GYRO, 
                             IDP_SAS_AGE_ACTIVITY, 
                             IDP_SAS_AGE_GRANULATION_CGBM)

    M_result = None
    if IDP_SAS_CONSISTENCY_FLAGS[0]:
        for M_name, M in zip(M_names, Ms):
            if M is None or len(M) <= 1:
                continue
            M_result = M_name
            break

    R_result = None
    if IDP_SAS_CONSISTENCY_FLAGS[1]:
        for R_name, R in zip(R_names, Rs):
            if R is None or len(R) <= 1:
                continue
            R_result = R_name
            break

    A_result = None
    if IDP_SAS_CONSISTENCY_FLAGS[2]:
        for A_name, A in zip(A_names, As):
            if A is None or len(A) <= 1:
                continue
            A_result = A_name
            break

    return [M_result, R_result, A_result]
