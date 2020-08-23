############################################################################# 80
'''
============================
fractoolbox Function Library
============================



'''

# Convert Fracture Data Format
# ============================
'''
Tools to convert between the various ways fracture geometry are described

--------------------
Data format examples 
--------------------
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

-------------------
The right-hand rule
-------------------
Place your right hand on the fracture plane with fingers pointed down
the direction of dip. Your thumb points to the strike azimuth which is
used under the right-hand rule convention. 
'''

def dip2strike(dipaz):
    '''Convert dip-dipazimuth data to strike using the right-hand rule
    #
    -----
    Input
    -----
    dipaz
        Azimuth of dip in degrees from north
    #
    ------
    Output
    ------
    return
        Strike azimuth in degrees from north 
        which adheres to the right hand rule 
    '''
    if dipaz < 90:
        strike=(dipaz-90)+360
    else:
        strike=dipaz-90
    return strike


def strike2dipaz(strike):
    '''Convert strike to dip azimuth using the right-hand rule convention
    #
    -----
    Input
    -----
    strike
        Strike azimuth in degrees from north 
        where the right-hand rule is observed
    #
    ------
    Output
    ------
    return
        Azimuth of the dip direction (0-360 degrees)
    '''
    if strike > 270:
        dipaz=(strike+90)-360
    else:
        dipaz=strike+90
    return dipaz

# Geometric Sample Bias: Isogenic Contours
# ========================================
'''

'''