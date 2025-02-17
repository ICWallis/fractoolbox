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

def simple_linear_Sv(maxdepth,obsdepth,density):
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



def Sv(depth, porosity, dry_rock_density, fluid_density):
    '''
    Calculate vertical stress in MPa for a given depth, porosity, dry rock density and fluid density.
        
    This function computes the vertical stress (Sv) at various depths using the bulk density profile 
    derived from input porosity, dry rock density, and fluid density. Vertical stress is calculated 
    through integration of gravitational force on the bulk density along depth.

    Parameters:
    - depth (array-like): Depth values in meters.
    - porosity (float or array-like): Porosity values (fraction of volume). If a single value is given, 
      it will be applied uniformly across all depths.
    - dry_rock_density (float or array-like): Density of dry rock in kg/m³. If a single value is given, 
      it will be applied uniformly across all depths.
    - fluid_density (float or array-like): Density of fluid in kg/m³. If a single value is given, 
      it will be applied uniformly across all depths.

    Returns:
    - vertical_stress_min_MPa (ndarray): Calculated vertical stress values at each input depth.
    
    Notes:
    - The function assumes gravitational acceleration to be 9.8 m/s².
    - The `integrate.cumtrapz` method is used for numerical integration.

    Suggestions:
    - Use the following units: depth in meters, density in kg/m³, porosity as a decimal percentage. 
    - Result Sv is returned in Pa. Multiply by 1e-6 to convert to MPa.

    '''
    # If dry_rock_density is a single value, convert to a list with the same length as depth
    if isinstance(dry_rock_density, (int, float)):
        dry_rock_density = [dry_rock_density] * len(depth)
    else:
        dry_rock_density = dry_rock_density
    
    # If porosity is a single value, convert to a list with the same length as depth
    if isinstance(porosity, (int, float)):
        porosity = [porosity] * len(depth)
    else:
        porosity = porosity

    # If fluid_density is a single value, convert to a list with the same length as depth
    if isinstance(fluid_density, (int, float)):
        fluid_density = [fluid_density] * len(depth)
    else:
        fluid_density = fluid_density
    
    
    # Calculate bulk density
    bulk_density = (1 - np.array(porosity)) * np.array(dry_rock_density) \
                             + np.array(porosity) * np.array(fluid_density)
    
    # Calculate gravitational force on bulk density and create arrays
    x = np.array(depth)
    y = np.array(bulk_density) * 9.8  # density array * gravity
    
    # Integrate arrays calculate vertical stress in Pa and convert to MPa
    vertical_stress = integrate.cumtrapz(y, x, initial=0)
    
    return vertical_stress


def estimate_shmin_cfc(Sv, Pp, mu):
    '''
    Use Coulomb frictional-failure theory for normal faulting to estimate Shmin.

    - Equation from Zoback (2010)
    
        Args:   Pp Pore pressure (MPa) 
                Sv Overburden stress (MPa)
                mu Coefficient of friction (unitless)
        
        Coefficient of friction should be around 0.6 in the deep reservior
        and around 0.3 in the clay cap and shallow aquifer zone 
            
        Returns: Estimate of Shmin (MPa)
    '''
    Shmin = ((Sv - Pp)/(((mu**2 + 1.)**0.5 + mu)**2)) + Pp
    return Shmin


def shmax_from_ESR(ESR_SHmax, Sv, Pp):
    '''calculate the magnitude of SHmax for a given SHmax Effective Stress Ratio, Pp and Sv.

        Parameters:
        -----------
        ESR_SHmax : float
            The Effective Stress Ratio for SHmax 
        Sv : float
            Vertical stress magnitude
        Pp : float
            Pore pressure magnitude

        Returns:
        --------
        SHmax : float
            The calculated magnitude of SHmax

        Formula:
        --------
        SHmax = ESR_SHmax * (Sv - Pp) + Pp

        Notes:
        ------
        - This calculation assumes a linear relationship between SHmax and effective stress as defined by the effective stress ratio.
        - Ensure consistent units are used for all inputs.
    '''
    SHmax = ESR_SHmax * (Sv - Pp) + Pp
    return SHmax

