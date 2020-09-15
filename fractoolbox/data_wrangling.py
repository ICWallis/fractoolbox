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
    '''Convert dip-dipazimuth data to strike using the right-hand rule

    Args:
        dipaz: Azimuth of dip in degrees from north
    
    Returns: 
        Strike azimuth (0-360 degrees) based on the right hand rule
 
    '''
    if dipaz < 90:
        strike=(dipaz-90)+360
    else:
        strike=dipaz-90
    return strike


def strike2dipaz(strike):
    '''Convert strike to dip azimuth using the right-hand rule convention

    Args:
        strike: Strike azimuth in degrees from north 
            where the right-hand rule is observed
 
    Returns
        Azimuth of the dip direction (0-360 degrees)
    
    '''
    if strike > 270:
        dipaz=(strike+90)-360
    else:
        dipaz=strike+90
    return dipaz

def xyzinterp(mDdat, mDsur, xsur, ysur, zsur):
    '''Interpolates xyz for a point on the well path
    
    Interpolation uses well survey data files and
    z can be replaced with any other attribute.
    Measured depth is the common value between the well survey data 
    and data that is receiving values (typically a dataframe of fractures)

    Args:
        xyzinterp expects input objects to be columns from Pandas dataframes
        Ensure all data has the same datum (typically rig floor for well data)

        mDdat: Meters along the well path of the data
        mDsur: Meters along the well path from the survey file
        xsur: Easting location from the survey file in whatever co-ordinate system 
        ysur: Northing location from the survey file in whatever co-ordinate system
        zsur: Vertical depth in total vertical depth (will not work with elevation data)
            or any other columns in the survey dataframe (e.g. plunge or azimuth)

    Returns:
        A dataframe with the mD, x, y, and z at the depth MD
        where z can be any data type
    
    Usage example:
        # append well plunge to fracture dataframe
        # dffracture: a pandas dataframe with fractures
        # dfsurvey: a pandas dataframe with well survey data 
    
        mDdat = dffracture['depth_mMDRF']
        mDsur = dfsurvey['depth_mMDRF']
        xsur = dfsurvey['easting_m']
        ysur = dfsurvey['northing_m']
        zsur = dfsurvey['plunge']

        dfxyz = ftb.xyzinterp(mDdat, mDsur, xsur, ysur, zsur) 
        dfxyz.columns = ['Dpth_mMDRT','Easting','Northing','wellplunge']
    
        # combine calculated values into fracture dataframe
        dffracture = pd.concat([dffracture,dfxyz1], axis=1, join='inner')
    
    '''
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


def linear_interpolate_2dp(depth, datadepth, data):
    '''
    For given depths, interpolate an array of values (to 2 decimal places) from and existing depth/data array
    
    Example of usage with Pandas
    depth = dffracture['DEPT']
    datadepth = dfsurvey['Depth_mMD']
    data = dfsurvey['plunge']
    dffracture['WellPlunge'] = linear_interpolate(depth, datadepth, data)
    '''
    values = np.around((np.interp(depth, datadepth, data)),2)
    return values


