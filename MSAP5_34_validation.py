#!usr/bin/env python 
#### PLATO MSAP5-34 Validation
#### Author: Earl Patrick Bellinger ( bellinger@phys.au.dk ) 
#### Max Planck Institute for Astrophysics, Garching, Germany 
#### Stellar Astrophysics Centre, Aarhus, Denmark 

import numpy as np
import pandas as pd
from scipy import stats
from MSAP5_32_selection import selection, M_names, R_names, A_names

models = pd.read_csv('MIST-models.csv')

def validation(IDP_MASS_SEISMIC=None, 
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
    Validates the selected combination of mass, radius, and age against a set of theoretical models.
    
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
    
    M_select, R_select, A_select = selection(*Ms, *Rs, *As)
    
    # compute means and stds 
    M_result = None
    if M_select is not None:
        M = Ms[M_names.index(M_select)]
        M_result = (M.mean(), M.std(), M_select)
    
    R_result = None
    if R_select is not None:
        R = Rs[R_names.index(R_select)]
        R_result = (R.mean(), R.std(), R_select)
    
    A_result = None
    if A_select is not None:
        A = As[A_names.index(A_select)]
        A_result = (A.mean(), A.std(), A_select)
    
    results = [M_result, R_result, A_result]
    fail = [None, None, None]
    
    # zero or only one measurement, no validation needed 
    if M_result is None and R_result is None or \
       M_result is None and A_result is None or \
       R_result is None and A_result is None:
        return results
    
    # validation: all three measurements 
    elif None not in results:
        chi2 = (((models['mass']   - M_result[0])/M_result[1])**2 + \
                ((models['radius'] - R_result[0])/R_result[1])**2 + \
                ((models['age']    - A_result[0])/A_result[1])**2)/3
        chi2min = chi2.min()
        print(stats.chi2.pdf(chi2min, 2))
        if stats.chi2.pdf(chi2min, 2) < 0.01: # reject the null hypothesis
            return fail
    
    # validation: two measurements 
    elif M_result is not None and R_result is not None:
        chi2 = (((models['mass']   - M_result[0])/M_result[1])**2 + \
                ((models['radius'] - R_result[0])/R_result[1])**2)/2
        chi2min = chi2.min()
        if stats.chi2.pdf(chi2min, 1) < 0.01: # reject the null hypothesis
            return fail
    elif M_result is not None and A_result is not None:
        chi2 = (((models['mass']   - M_result[0])/M_result[1])**2 + \
                ((models['age']    - A_result[0])/A_result[1])**2)/2
        chi2min = chi2.min()
        if stats.chi2.pdf(chi2min, 1) < 0.01: # reject the null hypothesis
            return fail
    elif R_result is not None and A_result is not None:
        chi2 = (((models['radius'] - R_result[0])/R_result[1])**2 + \
                ((models['age']    - A_result[0])/A_result[1])**2)/2
        chi2min = chi2.min()
        if stats.chi2.pdf(chi2min, 1) < 0.01: # reject the null hypothesis
            return fail
    
    return results
