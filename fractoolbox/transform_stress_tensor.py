

import numpy as np
import math
import pandas as pd



def Rs(alpha,beta,gamma):
    '''Generates an array that's used to transform (rotate) the stress tensor into a geographic coordinate system 
    
    Geographic coordinates are X North, Y East, and Z Down.
    Input Euler angles alpha, beta, gamma in degrees
    Defining the stress field in Euler angles...
    If S1 is vertical (normal faulting) then:
        alpha = the trend of SHmax - pi/2 (aka the azimuth in degrees minus 90 degrees)
        beta = the -ve trend of Sv (aka -90 for vertical stress)
        gamma = 0.
    If S1 is horizontal (strike slip or reverse faulting) then:
        alpha = trend of S1
        beta = -ve plunge of S1
        gamma = rake of S2
    Output is an array.
    Function is called by fSg that makes the matrix multiplication to do the transformation
    Method from appendix in Peska and Zoback (1995)
    '''
    alpha = math.radians(alpha)
    beta = math.radians(beta)
    gamma = math.radians(gamma)
    Rs = np.array([
                [np.cos(alpha)*np.cos(beta), 
                 np.sin(alpha)*np.cos(beta), 
                 -np.sin(beta)],
                [np.cos(alpha)*np.sin(beta)*np.sin(gamma)-np.sin(alpha)*np.cos(gamma), 
                 np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma), 
                 np.cos(beta)*np.sin(gamma)],
                [np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma),
                 np.sin(alpha)*np.sin(beta)*np.cos(gamma)-np.cos(alpha)*np.sin(gamma),
                 np.cos(beta)*np.cos(gamma)]
                ])
    return Rs 

def Rf(strike,dip):
    '''Generates an array for that's used to transform (rotate) the stress tensor from geographic to fracture/fault plane coordinates

    Input is strike and dip in degrees following the right hand rule
    (otherwise, if the fault dipped to the left when viewed along strike, then the dip would be a negitve number - not ideal)
    Returns a matrix that is used in the matrix multiplication that makes the conversion
    Function is called by the fSf that does the transformation 
    Method from pp 156-157 in Zoback (2010)'''
    strike = math.radians(strike)
    dip = math.radians(dip)
    Rf = np.array([
                    [np.cos(strike),np.sin(strike),0],
                    [np.sin(strike)*np.cos(dip), -np.cos(strike)*np.cos(dip), -np.sin(dip)],
                    [-np.sin(strike)*np.sin(dip),np.cos(strike)*np.sin(dip),-np.cos(dip)]
    ])
    return Rf

def rake(Sf):
    '''Boolen expression used to calculate the rake of a fracture
    
    Input the stress tensor in the coordinate system of the fracture
    Output is the rake of the fracture 
    Output is used in fRt to generate an array that transformations (rotates) the stress tensor into the the rake vector
    Function is called by fSr where the transformation into the rake vector occurs
    Contains optional print statements to show which statement is true
    Method from pp 156-157 in Zoback (2010)
    '''
    A = Sf[2,1]
    B = Sf[2,0]
    if A > 0.00001 and B > 0.00001:
        r = np.arctan(A/B) 
        #print('r is case 1')
    elif A > 0.00001 and B < 0.00001:
        r = np.arctan(A/B) 
        #print('r is case 2')
    elif A < 0.00001 and B >= 0.00001:
        r = math.radians(180)-np.arctan(A/-B)
        #print('r is case 3')
    elif A < 0.00001 and B < 0.00001:
        r = np.arctan(-A/-B)-math.radians(180)
        #print('r is case 4')
    return r

def Rt(rake):
    '''
    Generates an array for that's used to transform the stress tensor from fracture plane coordinates into the rake vector
    
    Input is the rake of the fracture/fault generated by frake
    Output is called by fSr used to transformation (rotate) the stress tensor into the frake vector where |Sr(3,1)| is the shear stress magnatude (tau)
    Method from pp 156-157 in Zoback (2010)
    '''
    Rt = np.array([
                [np.cos(rake),np.sin(rake),0],
                [-np.sin(rake),np.cos(rake),0],
                [0,0,1]
                ])
    return Rt



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
        
        Citation:   Zoback (2010) pp 156-157
    '''
    strike = math.radians(strike)
    dip = math.radians(dip)

    Rf = np.array([
        [np.cos(strike), np.sin(strike), 0],
        [np.sin(strike)*np.cos(dip), -np.cos(strike)*np.cos(dip), -np.sin(dip)],
        [-np.sin(strike)*np.sin(dip), np.cos(strike)*np.sin(dip), -np.cos(dip)]
    ])

    return Rf

def find_fracture_rake(Sf):
    '''Boolean expression used to calculate the rake of a fracture
    
    Input the stress tensor in the coordinate system of the fracture
    Output is the rake of the fracture 
    Output is used in fRt to generate an array that transformations (rotates) the stress tensor into the the rake vector
    Function is called by fSr where the transformation into the rake vector occurs
    Contains optional print statements to show which statement is true
    Method from pp 156-157 in Zoback (2010)
    '''
    A = Sf[2,1]
    B = Sf[2,0]
    if A > 0.00001 and B > 0.00001:
        r = np.arctan(A/B) 
        #print('r is case 1')
    elif A > 0.00001 and B < 0.00001:
        r = np.arctan(A/B) 
        #print('r is case 2')
    elif A < 0.00001 and B >= 0.00001:
        r = math.radians(180)-np.arctan(A/-B)
        #print('r is case 3')
    elif A < 0.00001 and B < 0.00001:
        r = np.arctan(-A/-B)-math.radians(180)
        #print('r is case 4')
    return r

def transform_tensor_from_fracture_plane_to_rake(rake):
    '''
    Generates an array for that's used to transform the stress tensor from fracture plane coordinates into the rake vector
    
    Input is the rake of the fracture/fault generated by frake
    Output is called by fSr used to transformation (rotate) the stress tensor into the frake vector 
    where |Sr(3,1)| is the shear stress magnitude (tau)
    Method from pp 156-157 in Zoback (2010)
    '''
    Rt = np.array([
                [np.cos(rake),np.sin(rake),0],
                [-np.sin(rake),np.cos(rake),0],
                [0,0,1]
                ])
    return Rt