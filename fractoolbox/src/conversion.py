# ================
# Conversion tools
# ================
'''
Convert between dip azimuth and strike azimuth using the right hand rule.

The Right-Hand Rule:

    Place your right hand on the fracture plane with fingers pointed down
    the direction of dip. Your thumb now points to the strike azimuth.

License
-------
The content of this repository is licensed under the Apache License, 
Version 2.0 (the "License"); you may use these files if you comply with 
the License, which includes attribution. You may obtain a copy of the 
License at http://www.apache.org/licenses/LICENSE-2.0

'''

def dip2strike(dipazimuth):
    """Convert dip azimuth to strike using the right-hand rule

    Args:
        dipaz (float): Azimuth of dip in degrees from north

    Returns:
        float: Strike azimuth (0-360 degrees) based on the right hand rule
    """
    if dipazimuth < 90:
        strike = (dipazimuth - 90) + 360
    else:
        strike = dipazimuth - 90
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