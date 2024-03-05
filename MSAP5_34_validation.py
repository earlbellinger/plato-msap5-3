#!usr/bin/env python 
#### PLATO MSAP5-34 Validation
#### Author: Earl Patrick Bellinger ( earl.bellinger@yale.edu ) 
#### Department of Astronomy, Yale University

import numpy as np
import pandas as pd
from scipy import stats
from MSAP5 import M_names, R_names, A_names

models = pd.read_csv('MIST-models.csv')

def validation(IDP_SAS_MASS_GRID_MIXED=None, 
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
               IDP_SAS_MASS_PRIORITY=None,
               IDP_SAS_RADIUS_PRIORITY=None,
               IDP_SAS_AGE_PRIORITY=None):
    """
    Validates the selected combination of mass, radius, and age against a set of theoretical models.
    
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

    # compute means and stds 
    M_result = None
    if IDP_SAS_MASS_PRIORITY is not None:
        M = Ms[M_names.index(IDP_SAS_MASS_PRIORITY)]
        mean, std = gumr(M.mean(), M.std())
        M_result = (mean, std, IDP_SAS_MASS_PRIORITY)
        #M_result = (M.mean(), M.std(), IDP_SAS_MASS_PRIORITY)
    
    R_result = None
    if IDP_SAS_RADIUS_PRIORITY is not None:
        R = Rs[R_names.index(IDP_SAS_RADIUS_PRIORITY)]
        mean, std = gumr(R.mean(), R.std())
        R_result = (mean, std, IDP_SAS_MASS_PRIORITY)
        #R_result = (R.mean(), R.std(), IDP_SAS_RADIUS_PRIORITY)
    
    A_result = None
    if IDP_SAS_AGE_PRIORITY is not None:
        A = As[A_names.index(IDP_SAS_AGE_PRIORITY)]
        mean, std = gumr(A.mean(), A.std())
        A_result = (mean, std, IDP_SAS_MASS_PRIORITY)
        #A_result = (A.mean(), A.std(), IDP_SAS_AGE_PRIORITY)
    
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
        pvalue = stats.chi2.pdf(chi2min, 2)
        if pvalue < 0.01: # reject the null hypothesis
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

def gumr(x_m, x_u):
    """
    Returns a correctly stated mean and standard deviation following 
    the Guide for the Expression of Uncertainty in Measurement. 
    """
    # Truncate x_u to two significant digits
    x_u_truncated = float(f"{x_u:.2g}")
    
    # Determine the number of digits after the decimal in x_u_truncated
    if '.' in str(x_u_truncated):
        decimal_places = len(str(x_u_truncated).split('.')[1])
    else:
        decimal_places = 0
    
    # Format x_m to match the number of digits after the decimal in x_u_truncated
    x_m_formatted = f"{x_m:.{decimal_places}f}"
    
    return float(x_m_formatted), x_u_truncated
