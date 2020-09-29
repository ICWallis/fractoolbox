# ============================
# fractoolbox Function Library
# ============================
'''
Tools to convert between the various ways fracture geometry are described

Structural Data Format Examples 
-------------------------------
The following examples all using the same fracture plane

    Borehole image log analysis:
    -   dip magnitude and dip azimuth
        Format example: dip = 30, dipaz = 350

    Structural geology:
    -   strike (0-360 degrees) and dip (0-90 degrees), using the right-hand rule
        Format example: 260/30

    -   strike (0-180 degrees), dip (0-90 degrees), and dip direction (East or West)
        Format example: 80/30/W

Refer to stereonet-basics.ipynb for visual examples, 
including how to describe and plot a line or rake

The Right-Hand Rule
-------------------
Place your right hand on the fracture plane with fingers pointed down
the direction of dip. Your thumb points to the strike azimuth which is
used under the right-hand rule convention. 

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
import mplstereonet
import pandas as pd

def dip2strike(dipaz):
    """Convert dip-dipazimuth data to strike using the right-hand rule

    Args:
        dipaz (float): Azimuth of dip in degrees from north

    Returns:
        float: Strike azimuth (0-360 degrees) based on the right hand rule
    """
    if dipaz < 90:
        strike = (dipaz - 90) + 360
    else:
        strike = dipaz - 90
    return strike


def strike2dipaz(strike):
    """Convert strike to dip azimuth using the right-hand rule convention

    Args:
        strike (float): Strike azimuth in degrees from north

    Returns:
        float: Azimuth of the dip direction (0-360 degrees)
    """
    if strike > 270:
        dipaz = (strike + 90) - 360
    else:
        dipaz = strike + 90
    return dipaz

def xyzinterp(mDdat, mDsur, xsur, ysur, zsur):
    """Interpolates xyz for a point on the well path
    
    Interpolation finds an xyz point on a line defined by an existing xy and measured depth
    z can be depth or any other float attribute (e.g., resistivity)
    Measured depth is the common value between the interpolation point and the line 
    Method was designed for adding xyz to a Pandas dataframe of fractures

    Args:
        xyzinterp expects input objects to be columns from a Pandas dataframe
        Ensure all depth data has the same datum (typically rig floor for well data)

        mDdat (float): Meters measured depth along the well path of the data point
        mDsur (float): Meters measured depth along the line (well path from a directional survey file)
        xsur (float): Easting location from the survey file in whatever co-ordinate system 
        ysur (float): Northing location from the survey file in whatever co-ordinate system
        zsur (float): Total vertical depth (or other data such as well plunge or azimuth) 
            assocated at the location defined by xsur and ysur. Will not work with z = elevation data, 
            so use 'datum - vertical depth' to find elevation of the point

    Returns:
        A dataframe with the mD (float), x (float), y (float), and z (float) at the point defined by mDdat
    
    Usage example:
        # append well plunge to fracture dataframe
        # dffracture: a Pandas dataframe with fractures
        # dfsurvey: a Pandas dataframe with well survey data 
    
        mDdat = dffracture['depth_mMDRF']
        mDsur = dfsurvey['depth_mMDRF']
        xsur = dfsurvey['easting_m']
        ysur = dfsurvey['northing_m']
        zsur = dfsurvey['plunge']

        dfxyz = ftb.xyzinterp(mDdat, mDsur, xsur, ysur, zsur) 
        dfxyz.columns = ['Dpth_mMDRT','Easting','Northing','wellplunge']
    
        # combine calculated values into fracture dataframe
        dffracture = pd.concat([dffracture,dfxyz1], axis=1, join='inner')
    
    """
    dfout = pd.DataFrame()
    dfout['MD'] = mDdat
    dfout['x'] = np.around(
        (np.interp
        (
        np.asarray(mDdat.tolist()), 
        np.asarray(mDsur.tolist()), 
        np.asarray(xsur.tolist())
        )
        ),2)
    dfout['y'] = np.around(
        (np.interp
        (
        np.asarray(mDdat.tolist()),
        np.asarray(mDsur.tolist()),
        np.asarray(ysur.tolist())
        )
        ),2)
    dfout['z'] = np.around(
        (np.interp
        (
        np.asarray(mDdat.tolist()), 
        np.asarray(mDsur.tolist()), 
        np.asarray(zsur.tolist())
        )
        ),2)
    return dfout


def linear_interpolate_2dp(depth, data, newdatadepth):
    """Interpolate data values (to 2 dp) for given depth(s) from and existing depth/data array

    Args:
        depth (float): List of depths for the 'data'
        data (float): List of data associated with 'depth'
        newdatadepth (float): List of depths that will be used to interpolate new data

    Returns:
        float: Data at the point defined by datadepth
    
    Usage example:
        Values are columns in a Pandas dataframe
        Ensure all depths are from the same dataum (typically rig floor)
    
        depth = dffracture['DEPT']
        datadepth = dfsurvey['Depth_mMD']
        data = dfsurvey['plunge']
        dffracture['WellPlunge'] = linear_interpolate(depth, datadepth, data)
    """
    values = np.around((np.interp(depth, datadepth, data)),2)
    return values