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

License 
-------
fractoolbox is distributed under an Apache 2.0 license
https://choosealicense.com/licenses/apache-2.0/
'''
import numpy as np
from scipy import integrate

def linear_Sv(maxdepth,obsdepth,density):
    """Magnitude of overburden stress [Sv in MPa] at a given observation depth

    Simple integration model that uses single average density and 
    returns Sv for the observation depth or list of depths.

    Args:
        maxdepth (float): The maximum depth of the stress model [m]
        obsdepth (float or list of floats): Depth(s) where Sv will be returned [m]
        density (float): average rock density [kg/m3] which is typically 2200 - 2800
        All args accept float or integer values

    Returns:
        Sv at obsdepth [MPa] as float or list of floats
    """
    depth_model = np.array([0,maxdepth])
    density_model = np.array([density,density])
    gravity = 9.8
    
    # trapezoid integration with unit conversion from Pa to MPa
    Sv_model = (integrate.cumtrapz(density_model * gravity, depth_model, initial=0)) * 1.e-6 
    
    # linear interpolation from the Sv model
    Sv_obsdepth = np.around((np.interp(obsdepth, depth_model, Sv_model)),2)
    
    return Sv_obsdepth
