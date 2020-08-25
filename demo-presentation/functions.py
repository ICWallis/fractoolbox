# Functions for use in these notebooks
import numpy as np

# Hoop stress
# -----------
def ftheta(n):
    '''generates a set of numbers in radians between 0 and 2pi (0 to 360 degrees equivlent)
    
    Used in many of the functions calculculting properties around the borehole wall.
    Input the number of numbers that will be generated at equal spacing between 0 and 360 degrees including the start and finish value.
    Output is a list of numbers in radians.
    By convention, theta is measured from the SHmax azumuth.'''
    theta=np.linspace(0,2*np.pi,n)
    return theta

def fSradialV(SHmax, Shmin, Pp, Pmud, R, r, theta):
    '''Calculates the magnatude of the radial stress around the wellbore in a vertical well.

    By convention is referred to as $\sigma_{rr}$.
    Note that when R = r, we are at the borehole wall and sigma_rr = deltaP
    Inputs:
    1. SHmax = magnatude of the maximum horozontal stress MPa (total stress not effective stress)
    2. Shmin = magnatude of the minimum horozontal stress MPa (total stress not effective stress)
    3. Pp = pore pressure in MPa
    4. Pmud = pressure inside the well in MPa (consider equivelent circulating density)
    5. R = wellbore radius
    6. r = depth of investigation as radial distance from the centre of the well
    7. theta = the azumuths around the wellbore signam_rr will be calculated for (refer to the ftheta function) 
    Output: a list of values for radial stress at the specified azumuths
    Reccomendation: plot these values on the y axis with theta on the x axis
    Fucntion was written by Evert using Kirsh (1898) as presented in Jager et al. (2007) and Zoback (2010)'''
    sigma_rr = (
                0.5*(SHmax+Shmin-2*Pp)
                *(1-(R/r)**2)
                +0.5*(SHmax-Shmin)
                *(1-4*(R/r)**2
                +3*(R/r)**4)
                *np.cos(2*theta)
                +(Pmud-Pp)
                *(R/r)**2
                )
    return sigma_rr

def fSeffhoopV(SHmax, Shmin, Pp, Pmud, sigma_Dt, R, r, theta):
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

def fSeffhoopmaxV (SHmax, Shmin, Pp, Pmud, sigma_Dt):
    '''Maximum effective hoop stres at the wellbore wall in a vertical well.

    By convetion is written as $\sigma_{\theta\theta}^{max}$
    Simplification of function fSeffhoopV that assumes R = r and a Shmin azumuth
    Input: 
    1. SHmax = magnatude of the maximum horozontal stress MPa (total stress not effective stress)
    2. Shmin = magnatude of the minimum horozontal stress MPa (total stress not effective stress)
    3. Pp = pore pressure in MPa
    4. Pmud = pressure inside the well in MPa (consider equivelent circulating density)
    5. sigma_Dt = magnatude of thermal stress MPa with deltaT = Twell - Tres (refer to the fsigma_Dt function)
    NOTE: Pa = MPa*1.e6
    deltaP (Pp - Pmud) has been checked and generates the expected reduction in stress if well pressure increases
    Output: a single value for maximum effective hoop stress
    Function witten by Irene using eq6.8 pp174 Zoback 2010'''
    sigma_ttmax = 3*SHmax - Shmin - 2*Pp - (Pp - Pmud) + sigma_Dt
    return sigma_ttmax


def fSeffhoopminV(SHmax, Shmin,  Pp, Pmud, sigma_Dt):
    '''Minimum effective hoop stres at the wellbore wall in a vertical well.

    By convetion is written as $\sigma_{\theta\theta}^{min}$
    Simplification of function fSeffhoopV that assumes R = r and a Shmax azumuth
    Input: 
    1. SHmax = magnatude of the maximum horozontal stress MPa (total stress not effective stress)
    2. Shmin = magnatude of the minimum horozontal stress MPa (total stress not effective stress)
    3. Pp = pore pressure in Pa
    4. Pmud = pressure inside the well in Pa (consider equivelent circulating density)
    5. sigma_Dt = magnatude of thermal stress Pa with deltaT = Twell - Tres (refer to the fsigma_Dt function)
    NOTE: Pa = MPa*1.e6
    deltaP (Pp - Pmud) has been checked and generates the expected reduction in stress if well pressure increases
    Output: a single value for maximum effective hoop stress
    Function witten by Irene using eq6.8 pp174 Zoback 2010'''
    sigma_ttmin = 3*Shmin - SHmax - 2*Pp - (Pp - Pmud) + sigma_Dt
    return sigma_ttmin

def fSfarfieldV(SHmax, Shmin, R, r, theta):
    '''Far field stresses (total stress, not effective) intersected by the vertical well.
    
    By convetion is written as $\sigma_{r\theta}$
    Input:
    1. SHmax = magnatude of the maximum horozontal stress MPa
    2. Shmin = magnatude of the minimum horozontal stress MPa
    3. R = wellbore radius
    4. r = depth of investigation as radial distance from the centre of the well
    when R = r, we are at the borehole wall
    5. theta = the azumuths around the wellbore signam_rr will be calculated for (refer to the ftheta function) 
    Output: a list of stress magnatudes at the angles of theta
    Reccomendation: plot these values on the y axis with theta on the x axis
    Fucntion was written by Evert using Kirsh (1898) as presented in Jager et al. (2007) and Zoback (2010)'''
    sigma_rt= (
                0.5*(SHmax-Shmin)
                *(1+2*(R/r)**2
                -3*(R/r)**4)
                *np.sin(2*theta)
                )
    return sigma_rt

def fSwbaxisV(SHmax, Sv, Shmin, nu, Pp, R, r, theta, sigma_Dt):
    '''Effective stress parallel to the wellbore axis in a vertical well.

    By convetion is written as $\sigma_{zz}$
    1. SHmax = magnatude of the maximum horozontal stress MPa (total stress not effective stress)
    2. Shmin = magnatude of the minimum horozontal stress MPa (total stress not effective stress)
    3. Sv = vertical stress (MPa)
    4. nu = Possions ratio (MPa)
    5. Pp = reservior pore pressure (MPa)
    6. R = wellbore radius (m)
    7. r = depth of investigation as radial distance from the centre of the well (m)
    Although, the unit of R and r doesnt matter, so long as they are the same.
    8. theta = azumuth around the wellbore wall
    9. sigma_Dt = magnatude of thermal stress Pa with deltaT = Twell - Tres (refer to the fsigma_Dt function) in MPa
    NOTE: Pa = MPa*1.e6
    Output: list of stress magnatudes, one for each theta value
    Written by Evert using Kirsh (1898) as written in Jager et al 2007 and Zoback 2010'''
    sigma_zz=(
        Sv-Pp
        -2*nu*(SHmax-Shmin)
        *(R/r)**2
        *np.cos(2*theta)
        +sigma_Dt)
    return sigma_zz

def fsigma_zz(sigma,nu,theta):
    '''some explnation'''
    sigma_one = sigma[0][0]        # sigma_11 MPa
    sigma_two = sigma[1][1]        # sigma_22 MPa
    sigma_three = sigma[2][2]      # simga_33 MPa
    sigma_tauA = sigma[0][1]       # sigma_12 MPa
    sigma_zz = (
            sigma_three 
            - 2*nu*(sigma_one-sigma_two) 
            * np.cos(2*theta) 
            - 4*nu*sigma_tauA*np.sin(2*theta)
            )
    return sigma_zz

def fsigma_tt(sigma,theta,deltaP):
    '''some explnation'''
    sigma_one = sigma[0][0]        # sigma_11 MPa
    sigma_two = sigma[1][1]        # sigma_22 MPa
    sigma_tauA = sigma[0][1]       # sigma_12 MPa
    sigma_tt = (
            sigma_one 
            + sigma_two 
            - 2*(sigma_one - sigma_two) 
            * np.cos(2*theta) 
            - 4*sigma_tauA*np.sin(2*theta) 
            - deltaP
            )
    return sigma_tt

def fsigma_tz(sigma,theta):
    '''some explnation'''
    sigma_tauB = sigma[1][2]       # sigma_23 MPa
    sigma_tauC = sigma[0][2]       # sigma_13 MPa
    tau_tz =   (
            2*(sigma_tauB*np.cos(theta) 
            - sigma_tauC*np.sin(theta))
            )
    return (tau_tz)

def fsigma_rr(deltaP):
    '''some explnation'''
    sigma_rr = deltaP
    return sigma_rr


# Thermal stress
# --------------
# look at updating this calculation using the Moos and Zoback paper

def fsigma_Dt(therex, K, nu, Tres, Twell):
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

# Vertical stress
# ---------------

def rohm(phi, rohr, rohf):
    '''mean rock density (rohm, kg/m3) that accounts for the effect porosity has on density

    porosity (phi, decimal %) 
    mineral density (rohr, kg/m3) 
    fluid density (rohf, kg/m3)
    Written by Irene using Jager et al 2007, pp 400'''
    rohm = (1-phi)*rohr+phi*rohf
    return rohm

def linSv(mdepth,obsdepth,dens):
    '''call the magnatued (MPa) of oberburden stress
    mdepth = bottom depth of the stress model in m
    obsdepth = the depths where Sv will be returned in m
    dens = rock density used in the model kg/m3
    this function assumes a single density with depth'''
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
