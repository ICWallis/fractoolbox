# Functions for use in these notebooks
import numpy as np
import pandas as pd
from scipy import integrate
from matplotlib import pyplot as plt


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
    
    Written by Irene using eq7.150  P204 Jager et al 2007
    
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

    Written by Evert using Kirsh (1898) as presented in Jager et al. (2007) and Zoback (2010)
    
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
    return sigma_tt


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

def minstress(S1,Pp,mu):
    '''Use the stress ratio and frictional faulting theroy to estimate the minimum stress
    
    Args:
        Pp: Pore pressure MPa 
        S1: Maximum stress MPa
        
    Returns:
        S3: Minimum stress    
    
    '''
    S3 = ((S1-Pp)/(((mu**2 + 1.)**0.5 + mu)**2))+Pp
    return S3

def maxstress(S3,Pp,mu):
    '''Using the stress ratio and frictional faulting theroy to estimate the maximum stress
    
    Args:
        Pp: Pore pressure MPa 
        S3: Minimum stress MPa
        
    Returns:
        S1: Maximum stress MPa
        
    '''
    S1 = ((S3-Pp)*(((mu**2 + 1.)**0.5 + mu)**2))+Pp
    return S1

def poly(Sv,Pp,mu,figname='StressPolygon'):
    '''Draws a stress polygon plot

    The stress polygon is the minimum and maximum horizontal stresses allowable based on the 
    Mohr-coloumb failure criterion. The edge of the polygon is where failure on an 
    optimally orientated fault or fracture will occur. 
    
    Args:
        Sv: Vertical stress [MPa]
        Pp: Pore pressure [MPa]
        mu: Coefficient of friction (0.1-0.6 with 0.5 a reasonable first estimate)

    Returns:
        A plot of stress polygons sometimes referred to as the Zoback-a-gram    
    '''
    minSh = minstress(Sv,Pp,mu)
    maxSh = maxstress(Sv,Pp,mu)
    minSH = minSh
    maxSH = maxSh
    ax = [minSh,minSh] # endpoints of the connecting lines
    ay = [minSh,Sv]
    bx = [minSh,Sv]
    by = [Sv,maxSH]
    cx = [Sv,maxSh]
    cy = [maxSH,maxSH]
    dx = [minSh,Sv]
    dy = [Sv,Sv]
    ex = [Sv,Sv]
    ey = [Sv,maxSH]
    fx = [minSh,maxSH]
    fy = [minSh,maxSH]
    f,ax1 = plt.subplots(1,1,figsize=(6,6))
    ax1.plot(ax,ay,color='k',alpha=0.5) # plots the connecting lines
    ax1.plot(bx,by,color='k',alpha=0.5)
    ax1.plot(cx,cy,color='k',alpha=0.5)
    ax1.plot(dx,dy,color='k',alpha=0.5)
    ax1.plot(ex,ey,color='k',alpha=0.5)
    ax1.plot(fx,fy,color='k',alpha=0.5)
    ax1.plot(minSh,minSH,'o',color='k')  # 1. Highly extensional $ S_hmin = S_Hmax << S_v $
    ax1.plot(Sv,Sv,'o',color='k')        # 2. Central point $ S_hmin = S_Hmax = S_v $
    ax1.plot(minSh,Sv,'o',color='k')     # 3. Transition between NF and SS $ S_hmin < S_Hmax = S_v $
    ax1.plot(Sv,maxSH,'o',color='k')     # 4. Transition between SS and RF $ S_hmin = S_v << S_Hmax$
    ax1.plot(maxSh,maxSH,'o',color='k')  # 5. Highy compresional $ S_hmin = S_Hmax >> S_v $
    ax1.text((maxSH-Sv)/4+Sv,(maxSH-Sv)/1.5+Sv, 'RF', fontsize=10)
    ax1.text((Sv-minSh)/2+minSh,(maxSH-Sv)/8+Sv, 'SS', fontsize=10)
    ax1.text((Sv-minSh)/4+minSh,(Sv-minSh)/1.5+minSh, 'NF', fontsize=10)
    plt.xlim(0,maxSh+20)
    plt.ylim(0,maxSh+20)
    plt.grid(linestyle='--')
    plt.xlabel('$S_{hmin}$ [MPa]')
    plt.ylabel('$S_{Hmax}$ [MPa]')
    plt.savefig(figname + '.png', dpi=300)
    plt.show()