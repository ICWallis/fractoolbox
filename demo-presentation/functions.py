# Functions for use in these notebooks
import numpy as np

def thermal_stress(therex, K, nu, Tres, Twell):
    '''Thermally induced stess [MPa] assuming a steady state has been reached.

    A convention of - as tensile and + as compressive has been used
    This convention means we use Twell-Tres here and that hoop stress 
    calculations above add Sigma_Dt.
    Inputs:
    therex = Coefficent of thermal expansion [typically 1.e-5 per kelvin]
    K = bulk modulus [typcilly 1.e10]
    nu = Possions ratio [typically 0.25]
    note that the elastic moduli need to be internally consistant
    Tres = reservior temp in kelvin
    Twell = well temp in Kelvin [typcally ~40degC in a geothermal well]
    celcius + 273.15 = Kelvin
    Output: a value of thermally induced stress.
    Written by Irene using eq7.150  P204 Jager et al 2007'''
    sigma_Dt = (
        (
            3*therex*K*
            ((1-2*nu)/(1-nu))
            *(Twell - Tres)
        )
        /1.e6)
    return sigma_Dt


def theta(n):
    '''generates a set of numbers in radians between 0 and 2pi (0 to 360 degrees equivlent)
    
    Used in many of the functions calculculting properties around the borehole wall.
    Input the number of numbers that will be generated at equal spacing between 0 and 360 degrees including the start and finish value.
    Output is a list of numbers in radians.
    By convention, theta is measured from the SHmax azumuth.'''
    theta=np.linspace(0,2*np.pi,n)
    return theta

def effhoopstress(SHmax, Shmin, Pp, Pmud, sigma_Dt, R, r, theta):
    '''Calculates the magnatude of the effective hoop stess around the wellbore in a vertical well.
    #
    By convention is referred to as $\sigma_{\theta\theta}$
    Note that when R = r, we are at the borehole wall
    #
    =============================
    Inputs
    =============================
    1. SHmax = magnatude of the maximum horozontal stress MPa (total stress not effective stress)
    2. Shmin = magnatude of the minimum horozontal stress MPa (total stress not effective stress)
    3. Pp = pore pressure in MPa
    4. Pmud = pressure inside the well in MPa (consider equivelent circulating density)
    5. sigma_Dt = magnatude of thermal stress MPa (refer to the fsigma_Dt function where deltaT = Twell - Tres)
    Becasue tension is here conceptuailsed as a -ve value (note how deltaT is calculated), we add sigma_Dt
    rather than subtracting as the equation appears in Zoback after Kirsh.
    6. R = wellbore radius
    7. r = depth of investigation as radial distance from the centre of the well
    8. theta = the azumuths around the wellbore signam_rr will be calculated for (refer to the ftheta function) 
    #
    =============================
    Output
    ============================= 
    a list of values for reffective hoop stess at the azumuths specified by theta
    #
    =============================
    Use Reccomendation 
    =============================
    Plot these values on the y axis with theta on the x axis
    =============================
    Dev History and Citation
    =============================
    Fucntion was written by Evert using Kirsh (1898) as presented in Jager et al. (2007) and Zoback (2010)'''
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
    return sigma_tt
