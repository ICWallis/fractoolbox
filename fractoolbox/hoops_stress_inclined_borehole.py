# Series of functions coded from Peska and Zoback (1995), with corrections from Zoback 2010,
# to model the maximum horizontal stress magnitude and orientation on an inclined borehole. 

import numpy as np
from matplotlib import pyplot as plt
import math

def normalise_by_value(input_list, value):
    result = []
    for n in input_list:
        result.append(n/value) 
    return result

def radians_to_degrees(radians):
    '''
    Convert list of radians to list of degrees

        Args: radians (list or 1D array) angle in radians

        Returns: list of given angles in degrees
    '''
    result = []
    for r in radians:
        degrees = math.degrees(r)
        result.append(degrees)
    return result


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
        [-np.cos(delta) * np.cos(phi), 
        -np.sin(delta) * np.cos(phi), 
        np.sin(phi)],
        [np.sin(delta), 
        -np.cos(delta),
        0],
        [np.cos(delta) * np.sin(phi),
        np.sin(delta) * np.sin(phi),
        np.cos(phi)]
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

        Args:       theta_degrees (list, float):  x axis value, azimuth around the borehole

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

    THis function has the assumptions that make S1 always vertical

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