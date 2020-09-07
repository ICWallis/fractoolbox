# ====================
# Stress Models
# ====================
'''In progress

'''
import numpy as np
import math
import pandas as pd
from scipy import integrate
from matplotlib import pyplot as plt

def linSv(mdepth,obsdepth,dens):
    '''Magnitude of overburden stress [Sv in MPa] at a given observation depth

    Integrates a single density with depth and then returns the a value
    from the curve at a desired depth of observation

    Args:
        mdepth = bottom depth of the stress model in m
        obsdepth = the depths where Sv will be returned in m
        dens = rock density used in the model kg/m3
        this function assumes a single density with depth

    Returns:
        Sv: Vertical stress     
    '''
    df=pd.DataFrame()                       # make a dataframe for the model
    df['depth'] = np.linspace(0,mdepth,100) # top and base depth of model
    df['dens'] = np.linspace(dens,dens,100)
    d = df['dens']                          # trapezoid integration
    g = 9.8            # gravity
    z = df['depth']
    x = z
    y = d*g   
    y_int = integrate.cumtrapz(y, x, initial=0)
    y_int_MPa_cons = y_int*1.e-6
    df['SvMPa'] = y_int_MPa_cons            # add Sv model dataframe for use later
    mDdata = obsdepth                       # grab the desired observation depth value
    mDsurvey = np.asarray(df['depth'].tolist())
    xsurvey = np.asarray(df['SvMPa'].tolist())
    Sv = np.around((np.interp(mDdata, mDsurvey, xsurvey)),2)
    return Sv
