# ============================
# fractoolbox Function Library
# ============================
'''
Python tools for structural geology and borehole image analysis that includes 
data handling, frequency and geometric analysis, and reservoir geomechanics.

Content
-------
Library content by section (and status)
-   Data Handling (in progress)
-   Geometric Sample Bias: Isogenic Contours (in progress)
-   Geomechanical Models (to come)
-   3DMohr Plot Analysis (to come)

Contributions
-------------
Initiated by Irene Wallis https://github.com/ICWallis/fractoolbox
as part of Doctoral Research at the University of Auckland that is 
supervised by David Dempsey https://github.com/ddempsey and 
Julie (JR) Rowland, with math/code contributions from Evert Durán 
https://github.com/edur409.

Citations
---------
Barton, C. A., Zoback, M. D., and Moos, D., 1995, Fluid flow along 
    potentially active faults in crystalline rock: Geology, v. 23, 
    p. 683-686.

Peška, P., and Zoback, M. D., 1995, Compressive and tensile failure 
    of inclined well bores and determination of in situ stress and 
    rock strength: Journal of Geophysical Research: Solid Earth, 
    v. 100, no. B7, p. 12791-12811.

Priest, S., 1993, Discontinuity Analysis for Rock Engineering, 
    Netherlands, Springer.

Terzaghi, R. D., 1965, Sources of error in joint surveys: Geotechnique, 
    v. 15, no. 3, p. 287-304.

Wallis, I.C., Rowland, J. V. and Dempsey, D. E., Allan, G., Sidik, R., 
    Martikno, R., McLean, K., Sihotang, M., Azis, M. and Baroek, M. 
    2020 (submitted) Approaches to imaging feedzone diversity with 
    case studies from Sumatra, Indonesia, and the Taupō Volcanic Zone, 
    New Zealand. New Zealand Geothermal Workshop: Waitangi, New Zealand.

Zoback, M. D., 2010, Reservoir Geomechanics, Cambridge University Presss

Licence 
-------
fractoolbox is distributed under an Apache 2.0 licence

https://choosealicense.com/licenses/apache-2.0/

'''
import numpy as np
import mplstereonet
import pandas as pd

# =============
# Data Handling
# =============
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
'''

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
        dfxyz.columns = ['Dpth_mMDRT','Northing','Easting','wellplunge']
    
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




