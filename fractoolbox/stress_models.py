# ========================
# Stress Tensor Estimation
# ========================
'''


Contributions
-------------
fractoolbox was initiated by Irene Wallis https://github.com/ICWallis/fractoolbox
as part of Doctoral Research at the University of Auckland that is 
supervised by David Dempsey https://github.com/ddempsey and 
Julie (JR) Rowland, with math/code contributions from Evert Durán 
https://github.com/edur409.

Licence 
-------
fractoolbox is distributed under an Apache 2.0 licence
https://choosealicense.com/licenses/apache-2.0/
'''
import numpy as np
from scipy import integrate

def linSv_3(maxdepth,obsdepth,density):
    '''Magnitude of overburden stress [Sv in MPa] at a given observation depth

    Simple intergration model using single density that returns a
    vertical stress vaule for the observation depth which is passed in.

    Args:
        maxdepth: the maximum depth of the stress model [m]
        obsdepth = depth(s) where Sv will be returned [m]
        obsdepth can be a single value, a list or a Pandas dataframe coloumn
        density = average rock density [kg/m3] which is typically 2200 - 2800

    Returns:
        Sv_obsdepth: Vertical stress aka overburden [MPa]     
    '''
    depth_model = np.array([0,maxdepth])
    density_model = np.array([density,density])
    gravity = 9.8
    # trapezoid integration with unit conversion from Pa to MPa
    Sv_model = (integrate.cumtrapz(density_model * gravity, depth_model, initial=0)) * 1.e-6 
    Sv_obsdepth = np.around((np.interp(obsdepth, depth_model, Sv_model)),2)
    return Sv_obsdepth
