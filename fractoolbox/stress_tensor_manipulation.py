# ====================================================================
# Stress Tensor Manipulation for Hoop Stress and Fracture Calculations
# ====================================================================
'''
Kirsh (1898) Equations for Stress Resolved onto a Circular Cavity

The Kirsh equations describe effective stresses around a wellbore in 
cylindrical coordinates for a vertical well and are used when the 
borehole axis = Sv. These can be used to explore the impact horizontal stress,
pore pressure and thermal stress have on the formation of drilling induced
borehole damage (i.e., borehole breakout and drilling induced tensile fractures).

Kirsch, 1898, Die theorie der elastizität und die bedürfnisse der festigkeitslehre: 
Zeitschrift des Vereines deutscher Ingenieure, v. 42, p. 797–807.

The functions below are from: 

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
fractoolbox is distributed under an Apache 2.0 license
https://choosealicense.com/licenses/apache-2.0/
'''

import numpy as np
from matplotlib import pyplot as plt
import math

#
# Make Stress Tensor
#


def make_stress_tensor(S1,S2,S3):
    '''Make an initial stress tensor

        Args:   S1  (float) Greatest principle stress, sigma 1
                S2  (float) Intermediate principle stress, sigma 2
                S3  (float) Least principle stress, sigma 3

        Returns: (2D array) an initial stress tensor

        NOTE function previously called fSs
    '''
    Ss = np.array([
        [S1,0,0],
        [0,S2,0],
        [0,0,S3]
    ])
    return Ss


def make_effective_stress_tensor(S1,S2,S3,Pp):
    '''Make an initial effective stress tensor

        Args:   S1  (float) Greatest principle stress, sigma 1
                S2  (float) Intermediate principle stress, sigma 2
                S3  (float) Least principle stress, sigma 3
                Pp  (float) Pore pressure

        Returns: (2D array) an initial effective stress tensor

        Notes:  Stresses and pressure are unit agnostic, but MPa is the ideal

                Pore pressure is multiplied with an identity matrix because
                effective stress only applies to volume change and, therefore, 
                not the shear stresses
        
        NOTE function previously called fSsEf

        Citation: Peska and Zoback (1995)
'''
    stress_tensor = np.array([
        [S1,0,0],
        [0,S2,0],
        [0,0,S3]
    ])
    # make identity matrix
    delta_ij = np.eye(3)
    # subtract pore pressure from normal stresses
    result = stress_tensor - (delta_ij * Pp)
    return result


#
# Transform Stress Tensor
#


def geographic_rotation_array(alpha,beta,gamma):
    '''Returns array used to transform an initial stress tensor into the geographic coordinate system 
    
        Args:   alpha (float) Euler angle alpha in degrees (refer to notes below)
                beta (float) Euler angle alpha in degrees (refer to notes below)
                gamma (float) Euler angle alpha in degrees (refer to notes below)

        Returns:    Array used for transformation
                    Referred to as Rs by Peska/Zoback 

        Notes:  Function called by transform_from_initial_to_geographic (fSg) which does the 
                matrix multiplication to execute the transformation

                Geographic coordinates are X North, Y East, and Z Down.
                
                Defining the stress field in Euler angles:
                    
                    If S1 is vertical (normal faulting) then:
                        alpha = the trend of SHmax - pi/2 (aka the azimuth in degrees minus 90 degrees)
                        beta = the -ve trend of Sv (aka -90 for vertical stress)
                        gamma = 0.
                    
                    If S1 is horizontal (strike slip or reverse faulting) then:
                        alpha = trend of S1
                        beta = -ve plunge of S1
                        gamma = rake of S2
                
                NOTE Function used to be called fRs

        Citation:   Peska and Zoback (1995)
    '''
    alpha = math.radians(alpha)
    beta = math.radians(beta)
    gamma = math.radians(gamma)

    Rs = np.array([
        [np.cos(alpha) * np.cos(beta), 
         np.sin(alpha) * np.cos(beta), 
         -np.sin(beta) ],
        [np.cos(alpha) * np.sin(beta) * np.sin(gamma) - np.sin(alpha) * np.cos(gamma), 
         np.sin(alpha) * np.sin(beta) * np.sin(gamma) + np.cos(alpha) * np.cos(gamma), 
         np.cos(beta) * np.sin(gamma) ],
        [np.cos(alpha) * np.sin(beta) * np.cos(gamma) + np.sin(alpha) * np.sin(gamma),
         np.sin(alpha) * np.sin(beta) * np.cos(gamma) - np.cos(alpha) * np.sin(gamma),
         np.cos(beta) * np.cos(gamma) ]
    ])
    return Rs 


def fracture_rotation_array(strike,dip):
    '''Generate array used to transform the stress tensor from geographic to fracture/fault plane coordinates

        Args:   strike (float)  Fracture strike in degrees
                dip (float)     Fracture dip in degrees

        Returns:    2D array used to do transformation
                    Referred to as Rf by Zoback (2010)

        Notes:  Strike and dip must follow the right hand rule

                If strike and dip doesn't follow the right hand rule and 
                the fault dipped to the left when viewed along strike, 
                then the dip would be a negative number (not ideal for handling here)

                Function is called by the fSf that does the transformation 

                NOTE function used to be called fRf
        
        Citation: Zoback (2010) pp 156-157
    '''
    strike = math.radians(strike)
    dip = math.radians(dip)

    Rf = np.array([
        [np.cos(strike), np.sin(strike), 0],
        [np.sin(strike)*np.cos(dip), -np.cos(strike)*np.cos(dip), -np.sin(dip)],
        [-np.sin(strike)*np.sin(dip), np.cos(strike)*np.sin(dip), -np.cos(dip)]
    ])

    return Rf


def borehole_rotation_array(well_azimuth,well_inclination):
    '''Returns array used to transform the initial stress tensor from geographic to borehole coordinates
    
        Args:   well_azimuth (float) Angle in degrees clockwise from north
                                Refereed to as 'delta' in Peska/Zoback
                well_inclination (float) Angle in degrees upward from vertical 
                                Otherwise referred to as drillers dip
                                Referred to as 'phi' in Peska/Zoback
    
        Returns:    2D Array used for transformation
                    Referred to as Rb in Peska/Zoback
    
        Notes:      This function is called by the fSb function which 
                    does the matrix multiplication to execute the stress
                    tensor transformation 

                    NOTE function used to be called fRb

        Citation:
                    Peska and Zoback (1995)
    '''
    delta = math.radians(well_azimuth)
    phi = math.radians(well_inclination)

    Rb = np.array([
        [-np.cos(delta) * np.cos(phi), -np.sin(delta) * np.cos(phi), np.sin(phi)],
        [np.sin(delta), -np.cos(delta), 0],
        [np.cos(delta) * np.sin(phi), np.sin(delta) * np.sin(phi), np.cos(phi)]
    ])

    return Rb


def transform_from_initial_to_borehole(initial_stress_tensor, alpha, beta, gamma, well_azimuth, well_inclination):
    '''Transform initial stress tensor in into borehole coordinates

    Transformation includes taking the intermediate step of rotating into 
    geographic coordinates using the Euler angles before rotating into
    borehole coordinates using the well azimuth and inclination.

        Args:   initial_stress_tensor (2D array) Either the principle stresses or effective stresses
                                Refer to functions 'make_stress_tensor' and 'make_effective_stress_tensor'
                alpha (float)   Euler angle alpha in degrees (refer to notes below)
                beta (float)    Euler angle alpha in degrees (refer to notes below)
                gamma (float)   Euler angle alpha in degrees (refer to notes below)
                well_azimuth (float)  Angle in degrees clockwise from north
                                Refereed to as 'delta' in Peska/Zoback
                well_inclination (float) Angle in degrees upward from vertical 
                                Otherwise referred to as drillers dip
                                Referred to as 'phi' in Peska/Zoback

        Returns:    Stress tensor with vectors rotated into borehole coordinates
                    Method referred to as Sb in Peska/Zoback

        Notes:  Functions called by this function 
                    geographic_rotation_array (referred to as Rs)
                    borehole_rotation_array (referred to as Rb)

                Geographic coordinates are X North, Y East, and Z Down.
                
                Defining the stress field in Euler angles:
                    
                    If S1 is vertical (normal faulting) then:
                        alpha = the trend of SHmax - pi/2 (aka the azumuth in degrees minus 90 degrees)
                        beta = the -ve trend of Sv (aka -90 for vertical stress)
                        gamma = 0.
                    
                    If S1 is horizontal (strike slip or reverse faulting) then:
                        alpha = trend of S1
                        beta = -ve plunge of S1
                        gamma = rake of S2
                
                NOTE This function used to be called fSb

        Citation:   Peska and Zoback (1995)
    '''
    Rs = geographic_rotation_array(alpha, beta, gamma)
    RsT = Rs.T
    Rb = borehole_rotation_array(well_azimuth, well_inclination)
    RbT = Rb.T
    result = Rb@RsT@initial_stress_tensor@Rs@RbT
    return result



#
# Hoop Stress
#
# Series of functions coded from Peska and Zoback (1995), with corrections from Zoback 2010,
# to model the maximum horizontal stress magnitude and orientation on an inclined borehole. 


def theta(n):
    """Generates a list of radian values between 0 and 2 pi (i.e., 0 to 360 degrees)
    
    Used in many of the functions calculating properties around the borehole wall.
    By convention, theta is measured from the SHmax azimuth.
    
    Args:
        n (float): The number of numbers that will be generated at equal spacing 
        between 0 and 360 degrees including the start and finish value.
    
    Returns:
        (float) A list of numbers in radians (typically referred to as theta)

    """
    theta = np.linspace(0, 2 * np.pi, n)
    return theta


def calculate_sigma_zz(theta, sigma_one, sigma_two, sigma_three, sigma_tauA, nu):
    '''
    Calculate 'vertical' stress on borehole (along the borehole axis)

        Args:   theta (list or 1D array)   Angle around borehole wall in radians
                sigma_one (float)   Postion 11 (normal stress) in the stress tensor rotated borehole coordinates
                sigma_two (float)   Postion 22 (normal stress) in the stress tensor rotated borehole coordinates
                sigma_three (float) Postion 33 (normal stress) in the stress tensor rotated borehole coordinates
                sigma_tauA (float)  Position 12 (shear stress) in the stress tensor rotated into borehole coordinates
                nu (float)          Possion's ratio
        
        Returns:    sigma_zz (list) Vertical stress on borehole at each given angle

        Notes:  Theta is measured from the low-side of the borehole in a deviated well
                In a vertical well, theta is conventionally measured from the azimuth or SHmax
                However, I am not sure if this convention for vertical wells applies here

                An initial stress tensor may be rotated into borehole coordinates 
                using the 'transform_from_initial_to_borehole' function

                sigma_one, sigma_two, sigma_3 and sigma_tauA are unit agnostic, but the convention is MPa

        Citation: Peska and Zoback (1995)
    '''
    result = []
    for t in theta:
        s_zz = (
            sigma_three - 2 * nu * (sigma_one - sigma_two) 
            * np.cos(2 * t) - 4 * nu * sigma_tauA * np.sin(2 * t)
            )
        result.append(s_zz)
    return result


def calculate_sigma_tt(theta, sigma_one, sigma_two, sigma_tauA, deltaP):
    '''
    Calculate hoop stress on borehole

        Args:   theta (list or 1D array)   Angle around borehole wall in radians
                sigma_one (float)   Postion 11 (normal stress) in the stress tensor rotated borehole coordinates
                sigma_two (float)   Postion 22 (normal stress) in the stress tensor rotated borehole coordinates
                sigma_tauA (float)  Position 12 (shear stress) in the stress tensor rotated into borehole coordinates
                deltaP (float)      Pressure difference between the borehole and reservoir
                                    (i.e., borehole fluid pressure minus reservoir pressure)

        Returns:    sigma_theta_theta (list) Hoop stress on borehole at each given angle

        Notes:  Theta is measured from the low-side of the borehole in a deviated well
                In a vertical well, theta is conventionally measured from the azimuth or SHmax
                However, I am not sure if this convention for vertical wells applies here

                An initial stress tensor may be rotated into borehole coordinates 
                using the 'transform_from_initial_to_borehole' function

                sigma_one, sigma_two and sigma_tauA are unit agnostic, but the convention is MPa

                NOTE Used in the theta_omega_printer function

        Citation: Peska and Zoback (1995)
    '''
    result = []
    for t in theta:
        s_tt = (
            sigma_one + sigma_two - 2 * (sigma_one - sigma_two) 
            * np.cos(2 * t) - 4 * sigma_tauA * np.sin(2 * t) - deltaP
            )
        result.append(s_tt)
    return result


def calculate_tau_tz(theta, sigma_tauB, sigma_tauC):
    '''
    Calculate shear stress on borehole

        Args:   theta (list or 1D array)   Angle around borehole wall in radians
                sigma_tauB (float)   Postion 22 (normal stress) in the stress tensor rotated borehole coordinates
                sigma_tauC (float)  Position 12 (shear stress) in the stress tensor rotated into borehole coordinates

        Returns:    tau_tz (list) Shear stress on borehole at each given angle

        Notes:  Theta is measured from the low-side of the borehole in a deviated well
                In a vertical well, theta is conventionally measured from the azimuth or SHmax
                However, I am not sure if this convention for vertical wells applies here

                An initial stress tensor may be rotated into borehole coordinates 
                using the 'transform_from_initial_to_borehole' function

                sigma_tauB and sigma_tauC are unit agnostic, but the convention is MPa

        Citation: Peska and Zoback (1995)
    '''
    result = []
    for t in theta:
        t_tz =   (2 * (sigma_tauB * np.cos(t) - sigma_tauC * np.sin(t)))
        result.append(t_tz)
    return result


def calculate_sigma_rr(theta, deltaP):
    '''
    Calculate hoop stress on borehole

        Args:   theta (list or 1D array)   Angle around borehole wall in radians
                deltaP (float)      Pressure difference between the borehole and reservoir
                                    (i.e., borehole fluid pressure minus reservoir pressure)

        Returns:    sigma_rr (list) Radial stress on borehole at each given angle

        Notes:  Theta is measured from the low-side of the borehole in a deviated well
                In a vertical well, theta is conventionally measured from the azimuth or SHmax
                However, I am not sure if this convention for vertical wells applies here

                An initial stress tensor may be rotated into borehole coordinates 
                using the 'transform_from_initial_to_borehole' function

                Pressure is unit agnostic but the convention is MPa

        Citation: Peska and Zoback (1995)
    '''
    result = []
    for t in theta:
        result.append(deltaP)
    return result


def calculate_omega_angle(sigma_zz, sigma_tt, tau_tz):
    '''
    Calculate angle omega between maximum effective principle stress at the borehole and the borehole axis

        Args:   sigma_zz (list or 1D array) 'Vertical' stress at borehole wall (along axis)
                sigma_tt (list or 1D array) Hoop stress at borehole wall
                tau_tz (list or 1D array) Shear stress at borehole wall

        Returns:    Omega (list) in degrees

        Notes:  All inputs were calculated in radians and this function converts to degrees

                Omega is the angle between the borehole axis and the maximum effective principle stress (sigma_tmax)
                that lies in a plane tangential to the borehole wall
                Omega is also the angle that a tensile fracture will be inclined relative to the borehole axis
                Omega + 90 is the angle between the borehole axis and minimum effective principle stress (sigma_tmin)

        Citation:   Peaka and Zoback (1995)
    '''
    result = []
    for zz, tt, tz in zip(sigma_zz, sigma_tt, tau_tz):
        o = (np.arctan(2 * tz / (zz - tt))) / 2
        o = math.degrees(o)
        result.append(o)
    return result


def calculate_sigma_tmax(sigma_zz, sigma_tt, tau_tz):
    '''
    Calculate maximum effective stress on plane tangential to the borehole wall (sigma_tmax)
    
        Args:   sigma_zz (list or 1D array) 'Vertical' stress at borehole wall (along axis)
                sigma_tt (list or 1D array) Hoop stress at borehole wall
                tau_tz (list or 1D array) Shear stress at borehole wall

        Returns: sigma_tmax (list) 

        Notes:  Function is unit agnostic, but the convention is MPa

        Citation: Peska and Zoback (1995)
    '''
    result = []
    for zz, tt, tz in zip(sigma_zz, sigma_tt, tau_tz):
        s_tmax = (0.5 * (zz + tt + np.sqrt((zz - tt)**2 + 4 * tz**2)))
        result.append(s_tmax)
    return result 


def calculate_sigma_tmin(sigma_zz, sigma_tt, tau_tz):
    '''
    Calculate minimum effective stress on plane tangential to the borehole wall (sigma_tmin)
    
        Args:   sigma_zz (list or 1D array) 'Vertical' stress at borehole wall (along axis)
                sigma_tt (list or 1D array) Hoop stress at borehole wall
                tau_tz (list or 1D array) Shear stress at borehole wall

        Returns: sigma_tmin (list) 

        Notes:  Function is unit agnostic, but the convention is MPa

        Citation: Peska and Zoback (1995)
    '''
    result = []
    for zz, tt, tz in zip(sigma_zz, sigma_tt, tau_tz):
        s_tmin = (0.5 * (zz + tt - np.sqrt((zz - tt)**2 + 4 * tz**2)))
        result.append(s_tmin)
    return result 


def thermal_stress(therex, K, nu, Tres, Twell):
    """Thermally induced stress [MPa] assuming a steady state has been reached
    
    In fractoolbox, tensile stress is -ve and compressive stress is +ve
    This convention means that Twell-Tres is used to find the deltaT 
    and that thermal stress is added to the hoop stress calculations 
    (eg the function calculate_effective_hoop_stress). This is the opposite convection to 
    to what is used in Zoback (2010) pp 174 eq 6.4.

    Args:   therex (float): Coefficient of thermal expansion, which is typically 1.e-5 per Kelvin
            K (float): Bulk modulus, which is typically 1.e10
            nu (float): Poisson's ratio, which is typically 0.25. 
            Tres (float): Reservoir temperature in Kelvin
            Twell (float): Internal well temperature in Kelvin (ie temp of injected fluid)

    Returns:    (float) Thermally induced stress, sigma_Dt
    
    Notes:  Function written by Irene using eq 7.150  P204 Jager et al (2007)
            Ensure that the elastic moduli (K & nu) are internally consistent

    """
    sigma_Dt = (
        (
            3 * therex * K 
            * ((1 - 2 * nu) / (1 - nu))
            * (Twell - Tres)
        )
        / 1.e6)
    return sigma_Dt


def calculate_effective_hoop_stress(SHmax, Shmin, Pp, Pmud, sigma_Dt, R, r, theta):
    """Calculates the magnitude of the effective hoop stress around a vertical borehole
    
    As a convention, effective hoop stress is referred to as $\sigma_{\theta\theta}$.
    
    Tension here conceptualized as a -ve value (note how deltaT is calculated in this function), 
    so we add sigma_Dt rather than subtracting as the equation appears in Zoback after Kirsh.
    
    Note that when R = r, we are at the borehole wall.
    
    Args:   SHmax (float): Magnitude of the maximum horizontal stress MPa (total stress not effective stress)
            Shmin (float): Magnitude of the minimum horizontal stress MPa (total stress not effective stress)
            Pp (float): Pore pressure in MPa
            Pmud (float): Pressure inside the well in MPa (consider equivalent circulating density)
            sigma_Dt (float): Magnitude of thermal stress MPa (refer to the thermal_stress function where deltaT = Twell - Tres)
            R (float): Wellbore radius
            r (float): Depth of investigation as radial distance from the centre of the well
            theta (float): the azimuths around the wellbore sigma_rr will be calculated for (refer to the theta function) 
    
    Returns:    (float) A list effective hoop stress at azimuths specified by theta (ie $sigma_{\tau\tau}$)

    Notes:  Function written by Evert using Kirsh (1898) as presented in Jager et al. (2007) and Zoback (2010)
            NOTE this is not the version of sigma_tt used in the theta_omega_printer function 
    """
    sigma_tt = (
                0.5 * (SHmax + Shmin - 2 * Pp)
                * (1 + (R/r)**2)
                - 0.5 * (SHmax - Shmin)
                * (1 + 3 * (R/r)**4)
                * np.cos(2 * theta)
                - (Pmud - Pp)
                * (R/r)**2
                + sigma_Dt
                )
    return sigma_tt


#
# Hoop stress plots
#


def normalise_by_value(input_list, value):
    result = []
    for n in input_list:
        result.append(n/value) 
    return result


def radians_to_degrees(radians):
    '''
    Convert list of radians to list of degrees

        Args:   radians (list or 1D array) angle in radians

        Returns: list of given angles in degrees
    '''
    result = []
    for r in radians:
        degrees = math.degrees(r)
        result.append(degrees)
    return result


def peska_plot(
    theta_degrees,
    sigma_zz_norm,
    sigma_tt_norm,
    tau_tz_norm,
    sigma_rr_norm,
    sigma_tmax_norm,
    sigma_tmin_norm,
    omega,
    theta_degrees_min_sigma_tmin,
    omega_min_sigma_tmin,
    ):
    '''
    Standardized plot for checking the model results

        Args:   theta_degrees (list, float):  x axis value, azimuth around the borehole

                # All of the following must have the same shape as theta_degrees
                # They are stresses projected into borehole coordinates
                sigma_zz_norm (list, float): vertical stress 
                sigma_tt_norm (list, float): hoop stress
                tau_tz_norm (list, float): shear stress
                sigma_rr_norm (list, float): radial stress
                sigma_tmax_norm (list, float): maximum normal stress, where the maximum point(s) on this curve defines the borehole breakout depending on the compressive rock strength
                sigma_tmin_norm (list, float): minimum normal stress, where the minimum point(s) on this curve defines the tensile fractures
                omega (list, float): Possible angles of the tensile fractures relative to the borehole axis
                
                # Single values plotted as vertical lines
                theta_degrees_min_sigma_tmin (float): Angle of the tensile fracture relative to the low side of the borehole
                omega_min_sigma_tmin (float): Angle of the tensile fracture relative to the borehole axis

        Returns:    A plot object

        Usage example:  plot = peska_plot(with args in here in the right order)
                        plt.savefig('ShmaxModelling-PeskaZoback-ModelApproach-ForwardModel.png', facecolor='w')
    '''
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(20,5))

    ax1.set_title('Effective stresses at borehole')

    ax1.plot(
        theta_degrees, 
        sigma_zz_norm,
        label = 'Vertical stress ($\sigma_{zz}$)'
    )

    ax1.plot(
        theta_degrees, 
        sigma_tt_norm,
        label = 'Hoop stress ($\sigma_{tt}$)'
    )

    ax1.plot(
        theta_degrees, 
        tau_tz_norm,
        label = 'Shear stress ($\Theta_{tz}$)'
    )

    ax1.plot(
        theta_degrees, 
        sigma_rr_norm,
        label = 'Radial stress ($\sigma_{rr} = \Delta P$)'
    )

    ax2.set_title('Effective stresses on the plane tangential to the borehole')

    ax2.plot(
        theta_degrees,
        sigma_tmax_norm,
        label = 'Max tangential stress ($\sigma_{tmax}$)'
    )

    ax2.plot(
        theta_degrees,
        sigma_tmin_norm,
        label = 'Min tangential stress ($\sigma_{tmin}$)'
    )

    ax3.set_title('Orientation of Sigma_tmax')

    ax3.plot(
        theta_degrees,
        omega,
        label = 'Orientation of $\sigma_{tmax}$ ($\omega$)'
    )

    for ax in [ax2, ax3]:
        ax.vlines(
            theta_degrees_min_sigma_tmin,
            -50,
            150,
            linestyle=':'
        )

        ax.hlines(
            omega_min_sigma_tmin,
            0,
            180,
            linestyle=':'
        )

    for ax in [ax1, ax2, ax3]:
        ax.legend()
        ax.set_xlabel('Borehole azimuth from low-side [degrees]')
        ax.set_xlim(0,180)
        ax.minorticks_on()
        ax.grid(which='major',color='k',linewidth=0.5,alpha=0.5)
        ax.grid(which='minor', color='k', linestyle=':', linewidth=0.5,alpha=0.5)

    for ax in [ax1, ax2]:
        ax.set_ylabel('Effective stress / Sigma 1')
        ax.set_ylim(-3, 10)

    ax3.set_ylim(-50, 150)
    ax3.set_ylabel('Azimuth of Omega in degrees')

    return plt


def theta_omega_printer(
        S1,
        S2s, # function designed so more than one can be passed, must be list
        S3,
        Pp, 
        SHmax_orientations, # function designed so more than one can be passed, must be lists
        beta, 
        gamma, 
        well_az, 
        well_inc,
        deltaP,
        nu,
        ):
    '''
    Function that pints the theta (distance from low side)
    and omega (tilt) of the tensile fracture

    This function assumes that S1 always vertical
    S2 cases are varied to generate SS/RF conditions
    '''
    # make empty lists to append thetas and omegas to
    thetas = []
    omegas = []

    for S2, SHmax_ort in zip(S2s, SHmax_orientations):
        alpha = SHmax_ort - 90
        # =======================================================================================
        # Forward model: Given a known S2 and alpha, what is the tensile fracture theta and omega
        # =======================================================================================
        # Implemented as a two step process that brute forces the result
        # Step 1: Makes a list with all of the possible options
        # Step 2: Find the minimum theta of that list and the index location of that minimum value
        # Step 3: Use that index value to call the omega 
        # Step 4: Plot the results as a sense check

        # -----------------------------------------------------
        # Step 1: Makes a list with all of the possible options
        # -----------------------------------------------------

        # Make effective stress tensor
        effective_stress_tensor = make_effective_stress_tensor(S1,S2,S3,Pp)

        # Rotate stress tensor into geographic and then borehole coordinates
        effective_stress_tensor_borehole = transform_from_initial_to_borehole(effective_stress_tensor, alpha, beta, gamma, well_az, well_inc)

        # Make objects from the stress tensor                     # sigma_ij
        sigma_one = effective_stress_tensor_borehole[0][0]        # sigma_11 MPa
        sigma_two = effective_stress_tensor_borehole[1][1]        # sigma_22 MPa
        sigma_three = effective_stress_tensor_borehole[2][2]      # sigma_33 MPa
        sigma_tauA = effective_stress_tensor_borehole[0][1]       # sigma_12 MPa
        sigma_tauB = effective_stress_tensor_borehole[1][2]       # sigma_23 MPa
        sigma_tauC = effective_stress_tensor_borehole[0][2]       # sigma_13 MPa

        # Theta a around the borehole for calculations (radians)
        theta = np.linspace(0, math.radians(180), 3600) 
        # Only half the borehole done because the other half is the same

        # Theta for plotting (degrees)
        theta_degrees = radians_to_degrees(theta)

        # Calculate stresses on borehole wall
        sigma_zz = calculate_sigma_zz(theta, sigma_one, sigma_two, sigma_three,  sigma_tauA, nu)
        sigma_tt = calculate_sigma_tt(theta, sigma_one, sigma_two, sigma_tauA, deltaP)
        tau_tz = calculate_tau_tz(theta, sigma_tauB, sigma_tauC)
        sigma_rr = calculate_sigma_rr(theta, deltaP)

        # Calculate stresses on the plan tangential to the borehole wall
        sigma_tmax = calculate_sigma_tmax(sigma_zz, sigma_tt, tau_tz)
        sigma_tmin = calculate_sigma_tmin(sigma_zz, sigma_tt, tau_tz)

        # Calculate the angle between Sigma_tmax and the borehole axis
        omega = calculate_omega_angle(sigma_zz, sigma_tt, tau_tz)

        # Normalize results to S1
        sigma_zz_norm = normalise_by_value(sigma_zz, sigma_one)
        sigma_tt_norm = normalise_by_value(sigma_tt, sigma_one)
        tau_tz_norm = normalise_by_value(tau_tz, sigma_one)
        sigma_rr_norm = normalise_by_value(sigma_rr, sigma_one)
        sigma_tmax_norm = normalise_by_value(sigma_tmax, sigma_one)
        sigma_tmin_norm = normalise_by_value(sigma_tmin, sigma_one)

        # ----------------------------------------------------------------------------------------
        # Step 2: Find the minimum theta of that list and the index location of that minimum value
        # ----------------------------------------------------------------------------------------

        # Find the location in the list of minimum sigma_tmin and 
        # use that index location to find the corresponding theta angle
        min_sigma_tmin = min(sigma_tmin)
        index_min_sigma_tmin = sigma_tmin.index(min_sigma_tmin)
        tensile_fracture_theta = theta_degrees[index_min_sigma_tmin]

        #tensile_fracture_theta 
        #print('Calculated tensile fracture theta = ', tensile_fracture_theta)
        thetas.append(tensile_fracture_theta)
        #print('Case study tensile fracture theta = ', tensile_fracture_theta_measured)
        # theta_degrees_min_sigma_tmin is the angle in degrees from the low side of the borehole
        # to the center of the tensile fracture. This will be observed data.

        # For those non-unique cases where there is more than one minimum
        # the method above will only return one of the mins
        # Is there a way to find the other min? Must be some kind of root finding problem.
        # Or just find a couple of the lowest values

        # ----------------------------------------------
        # Step 3: Use that index value to call the omega 
        # ----------------------------------------------

        # Use the index location of the minimum sigma_tmin to to find the angle omega 
        tensile_fracture_omega = omega[index_min_sigma_tmin]
        #print('Calculated tensile fracture Omega = ', tensile_fracture_omega)
        omegas.append(tensile_fracture_omega)
        #print('Case study tensile fracture Omega = ', tensile_fracture_omega_measured)

        # omega_min_sigma_tmin is the angle of the maximum effective stress on the plane tangential to the borehole wall
        # which is also the angle of the tensile fracture relative to the borehole axis (one of our knows)
        # omega + 90 is the minimum principle stress on the tangential plane, and the value often used by peska
    return thetas, omegas


#
# Fracture Stress Calculations
#


def find_fracture_rake(Sf):
    '''Boolean expression used to calculate the rake of a fracture
    
    Input the stress tensor in the coordinate system of the fracture.
    Output is the rake of the fracture. 
    
    Output is used in transform_tensor_from_fracture_plane_to_rake (Rt) to generate an array that 
    transformations (rotates) the stress tensor into the the rake vector.
    
    Function is called by fSr where the transformation into the rake vector occurs

    Contains optional print statements to show which statement is true
    Method from pp 156-157 in Zoback (2010)
    '''
    A = Sf[2,1]
    B = Sf[2,0]
    if A > 0.00001 and B > 0.00001:
        rake = np.arctan(A/B) 
        #print('r is case 1')
    elif A > 0.00001 and B < 0.00001:
        rake = np.arctan(A/B) 
        #print('r is case 2')
    elif A < 0.00001 and B >= 0.00001:
        rake = math.radians(180)-np.arctan(A/-B)
        #print('r is case 3')
    elif A < 0.00001 and B < 0.00001:
        rake = np.arctan(-A/-B)-math.radians(180)
        #print('r is case 4')
    return rake


def transform_tensor_from_fracture_plane_to_rake(rake):
    '''
    Generates an array for that's used to transform the stress tensor from fracture plane coordinates into the rake vector
    
    Input is the rake of the fracture/fault generated by the function find_fracture_rake.
    Output is called by fSr used to transformation (rotate) the stress tensor into the rake vector.

    where |Sr(3,1)| is the shear stress magnitude (tau)

    Method from pp 156-157 in Zoback (2010)
    '''
    Rt = np.array([
                [np.cos(rake),np.sin(rake),0],
                [-np.sin(rake),np.cos(rake),0],
                [0,0,1]
                ])
    return Rt


# Template for the use of the funcations above

def fracture_sn_tau(S1,S2,S3,Pp,Norm,alpha,beta,gamma,strike,dip):
    '''Calculate the shear (tau) and normal (Sn) stress on a fracture

    Normalised to either vertical stress or effective vertical stress

    Args:
        S1:
        S2:
        S3:
        Pp:
        Sv:
        alpha:
        beta:
        gamma:
        strike:
        dip:

        Recommendation: this can be efficiently done using a tuple 
    
    Returns: 
        Sn: Stress normal to the fracture plane [MPa]
        tau: Shear stress on the fracture plane [MPa]

    '''
    '''
    # create effective stress array
    Ss = np.array([                 
                    [S1-Pp,0,0],
                    [0,S2-Pp,0],
                    [0,0,S3-Pp]
                    ])
    #print('Ss: the effective stress tensor =','\n',Ss,'\n')

    # use the three Euiler angles to generate an array that is used 
    # to transform the effective stress array into geographic coordinates
    # x has been added to Rs to differentiate it to Rs above
    Rsx = Rs(alpha,beta,gamma)  
    #print('Rs: the stress coordinate system based on' + 
    #   'the inputted Euler angles =','\n',Rsx,'\n')

    # rotate the stress tensor into geographic coordinates
    Sg = Rsx.T@Ss@Rsx                 
    #print('Sg: the effective stress tensor now rotated into' +
    #   'geographic coordinates =','\n',Sg,'\n')

    # use the fracture strike an dip to generate an array that is used
    # to transform the stress field into fracture coordinates
    Rfx = Rf(strike,dip)        
    #print('Rf: the matrix that rotates the stress tensor from' + '
    #   'geographic coordinates into the fracture plane coordinates =','\n',Rf,'\n')

    # transform the stress field into the fracture coordinate system
    # x has been added to Rf to differentiate it to Rf above
    Sf = Rfx@Sg@Rfx.T                 
    #print('Sf: the effective stress tensor now rotated into the' +'
    #   fracture plane coordinate system =','\n',Sf,'\n')

    # take stress normal to the fault plane and normalise it to 
    # vertical stress so fractures from different depths can be plotted together
    Sn = Sf[2,2]/Norm                 
    #print('Sn: the effective stress magnitude normal to the' + '
    #   'fault plane (Sf bottom right) normalised to Sv =','\n',Sn,'\n')

    # calculate the rake of the fracture assuming only dip-slip
    # x has been added to rake to differentiate it from rake function above
    rakex = rake(Sf)            
    #print('the rake of the slip vector =',rake,'\n')

    # use the rake to generate an array to transform the stress field into the rake
    Rtx = Rt(rakex)
    #print('Array Rt that is used to transform the effective stress from the' + 
    #   'fault co-ordinate system into the rake coordinate system so we can' + 
    #   'grab the shear stress magnitude along the slip vector =','\n',Rt,'\n')

    # transform the stress field into the direction of the rake
    # x has been added to Rt to differentiate it from the function above
    Sr = Rtx@Sf@Rtx.T
    #print('Effective stress tensor along the slip vector (Sr) where the' + 
    #   'bottom left (as an absolute number) is the shear stress magnitude (tau) =','\n',Sr,'\n')

    # take the absolue number of shear stress in the direction of the rake 
    # and normalise to vertical stress for plotting
    tau = abs(Sr[2,0])/Norm
    #print('Shear stress on the plane (tau, which is bottom left of Sr array) =',tau,'\n')

    return (Sn,tau)
    '''