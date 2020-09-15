# =================================================================
# Kirsh (1898) Equations for Stress Resolved onto a Circular Cavity
# =================================================================
'''
The Kirsh equations describe effective stresses around a wellbore in 
cylindrical coordinates for a vertical well and are used when the 
borehole axis = Sv. These can be used to explore the impact horizontal stress,
pore pressure and thermal stress have on the formation of drilling induced
borehole damage (i.e., borehole breakout and drilling induced tensile fractures).

Kirsch, 1898, Die theorie der elastizität und die bedürfnisse der festigkeitslehre: 
Zeitschrift des Vereines deutscher Ingenieure, v. 42, p. 797–807.

The functions in this notebook are the form of the equations as presented 
Jaeger et al (2007) and Zoback (2010). 

Jaeger, J. C., 2007, Fundamentals of rock mechanics, Malden, MA, Malden, 
MA : Blackwell Pub. 2007.

Zoback, M. D., 2010, Reservoir Geomechanics, Cambridge University Press.

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


def thermal_stress(therex, K, nu, Tres, Twell):
    '''Thermally induced stess [MPa] assuming a steady state has been reached.

    A convention of - as tensile and + as compressive has been used
    This convention means we use Twell-Tres here and that hoop stress 
    calculations add Sigma_Dt (thermal stress).
    
    Args:
        therex: Coefficient of thermal expansion [typically 1.e-5 per kelvin]
        K: bulk modulus [typically 1.e10]
        nu: Possions ratio [typically 0.25]
        Ensure that the elastic moduli need to be internally consistent
        Tres: Reservoir temp in kelvin
        Twell: Well temp in Kelvin [typically ~40degC in a geothermal well]

    Returns: 
        sigma_Dt: Thermally induced stress
    
    Function written by Irene using eq7.150  P204 Jager et al 2007
    
    '''
    sigma_Dt = (
        (
            3*therex*K*
            ((1-2*nu)/(1-nu))
            *(Twell - Tres)
        )
        /1.e6)
    return sigma_Dt

def theta(n):
    '''generates a set of numbers in radians between 0 and 2pi (0 to 360 degrees equivalent)
    
    Used in many of the functions calculating properties around the borehole wall.
    By convention, theta is measured from the SHmax azimuth.
    
    Args:
        n: The number of numbers that will be generated at equal spacing 
        between 0 and 360 degrees including the start and finish value.
    
    Returns:
        theta: A list of numbers in radians

    '''
    theta=np.linspace(0,2*np.pi,n)
    return theta

def effhoopstress(SHmax, Shmin, Pp, Pmud, sigma_Dt, R, r, theta):
    '''Calculates the magnitude of the effective hoop stress around the wellbore in a vertical well.
    
    By convention is referred to as $\sigma_{\theta\theta}$
    Because tension is here conceptualised as a -ve value (note how deltaT is calculated in this function), 
    we add sigma_Dt rather than subtracting as the equation appears in Zoback after Kirsh.
    Note that when R = r, we are at the borehole wall
    
    Args:
        SHmax: Magnitude of the maximum horizontal stress MPa (total stress not effective stress)
        Shmin: Magnitude of the minimum horizontal stress MPa (total stress not effective stress)
        Pp: Pore pressure in MPa
        Pmud: Pressure inside the well in MPa (consider equivalent circulating density)
        sigma_Dt: Magnitude of thermal stress MPa (refer to the fsigma_Dt function where deltaT = Twell - Tres)
        R: Wellbore radius
        r: Depth of investigation as radial distance from the centre of the well
        theta: the azimuths around the wellbore signam_rr will be calculated for (refer to the ftheta function) 
    
    Returns:
        sigma_tt: A list effective hoop stress at azimuths specified by theta (sigma tau tau)

    Function written by Evert using Kirsh (1898) as presented in Jager et al. (2007) and Zoback (2010)
    
    '''
    sigma_tt = (
                0.5*(SHmax+Shmin-2*Pp)
                *(1+(R/r)**2)
                -0.5*(SHmax-Shmin)
                *(1+3*(R/r)**4)
                *np.cos(2*theta)
                -(Pmud-Pp)
                *(R/r)**2
                +sigma_Dt
                )
    return sigma_t