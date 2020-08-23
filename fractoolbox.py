############################################################################# 80
'''
============================
fractoolbox Function Library
============================



'''

# Convert Fracture Data Format
# ============================
'''


-------------------
The Right-Hand Rule
-------------------
Place your right hand on the fracture plane with fingers pointed down
the direction of dip. Your thumb points to the strike azumuth which is
used under the right-hand rule convention. 
'''

def dip2strike(dipaz):
    '''Convert dip-dipazimuth data to strike using the right-hand rule
    #
    -----
    Input
    -----
    dipaz
        Azumuth of dip in degrees from north
    #
    ------
    Output
    ------
    return
        Strike azumuth in degrees from north 
        which adhears to the right hand rule
    '''
    if dipaz < 90:
        strike=(dipaz-90)+360
    else:
        strike=dipaz-90
    return strike

def strike2dipaz(strike):
    '''Convert strike to dip azumuth using the right-hand rule convention
    #
    -----
    Input
    -----
    strike
        Strike azumuth in degrees from north 
        where the right-hand rule is observed
    #
    ------
    Output
    ------
    return
        Azumuth of the dip direction (0-360 degrees)
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