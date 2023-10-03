# ========================================
# Geometric Sample Bias: Isogenic Contours
# ========================================
'''
Seminal work by Terzaghi (1965) revealed a geometric bias is generated 
by sampling a three-dimensional fracture network with a line (well path or scan-line). 
Fracture planes perpendicular to the line are likely to be intersected while 
those parallel to the line are rarely intersected. This geometric sample bias 
generates a 'blind zone' where fractures near-parallel to the line are missing.

Terzaghi (1965) proposed a methodology quantifies the geometric sample 
bias using the acute angle (alpha) between the fracture plane and 
the line. Visualizing the blind zone (where sin(alpha) +/- 0.3) and 
contours of sample bias (isogenic contours) on a stereonet enables
us to visually evaluate the degree that geometric sample bias in 
affects a fracture dataset.  

A weighting may be applied based on the alpha angle which may correct
the sampled fracture population to something more reflective of actual
frequency. This kind of correction is common-place in modern image log
analysis and some form of the Terzaghi correction comes baked into 
most log analysis software. However, there are two key issues with these
corrections:

-   Correction may mislead interpretation by emphasizing solitary 
    fractures that are not part of some significant but under-sampled 
    population, especially where the weighting factor approaches 
    infinity near sin α = 0. Priest (1993) recommends resolving this 
    by using an upper limit of sin α = 0.1 when weighting.

-   Correction can only be applied to those fractures which were sampled
    and therefore does a poor job of correcting in the blind zone where
    fractures are rarely sampled. 

Functions are included below that (1) weight fracture populations to 
reduce the impact of geometric sample bias and (2) construct isogenic 
contours for stereonets so the effect of sample bias and the blind 
zone are visible to the interpreter. Refer to Wallis et al. (2020)
for examples of these methods as applied micro-resistivity image 
logs acquired in seven high-temperature geothermal wells.

Refer to geometric-sample-bias.ipynb for an illustration of the 
alpha angles, code examples using the functions, and plots that show
the impact of the blind zone on a range of well data. 

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

# The geometry in these unit vector functions are adapted from 
# Priest (1993) to make the isogenic contours method work. The
# difference is probably abstraction issue where Preist uses 
# the pole (plane normal) throughout his book and I default to 
# dip-azimuth and well-azimuth, as is typical in borehole 
# analysis. I may at some point refactor the isogenic contour 
# method and change these vector functions. 

def unitvectorx(az,pl):
    '''Calculate the unit vector for a pole or well path
    
    Args: 
        az: Fracture dip azimuth or well azimuth
        pl: Plunge
    
    Returns: The unit vector x component

    '''
    x = -1 * np.sin(np.deg2rad(az)) * np.cos(np.deg2rad(pl))
    return x


def unitvectory(az,pl):
    '''Calculate the unit vector for a pole or well path

    Args:
        az: Fracture dip azimuth or well azimuth
        pl: Plunge
    
    Returns: The  unit vector y component

    '''
    y = -1 * np.cos(np.deg2rad(az)) * np.cos(np.deg2rad(pl))
    return y


def unitvectorz(az,pl):
    '''Calculate the unit vector for a pole or well path
    
    Args:
        az: Fracture dip azimuth or well azimuth
        pl: Plunge
    
    Returns: The  unit vector z component

    '''
    z = 1 * np.sin(np.deg2rad(pl))
    return z


def isogeniccontour(wpl, waz, sin_al):
    '''
    Parametrically calculates isogenic contours based on well path

    Args:
        wpl: Well plunge (0-90 degrees) which is the inclination of the
            well path measured from a horizontal surface
            
            Suggestion: If the well deviation is relativity consistent, 
                calculate the mean deviation magnitude (drillers dip) and 
                subtract this angle from 90 to find the plunge. If there
                are large variations in well geometry, then fracture data
                may need to be split onto separate stereonets, each with a
                set of isogenic contours calculated using that section's 
                mean well deviation. 

        waz: Well azimuth (0-360 degrees) which is the map direction
            the well is deviated in

            Suggestion: If the well azimuth is relativity consistent,
                calculate the mean azimuth. If there are large variations
                in well geometry, then proceeded as per suggestion made 
                above for varying well plunge.

        sin_al: Sin(alpha) (0.01 - 1) which is the contour intervals
            that will ascribe on the stereonet the likelihood a 
            fracture will be sampled

            Suggestion: a float at 0.1 increments from 0.1 to 0.9

    Returns: 
        strike: Strike azimuth (0-360 degrees) based on the right hand rule
        
        dip: Dip magnitude (0-90 degrees)
    
    Function developed by David Dempsey and Irene Wallis
    '''
    # Convert well azimuth/plunge to a vector p
    p = np.array([f(waz,wpl) for f in [unitvectorx, unitvectory, unitvectorz]])
    # make unit vector n1 that is perpendicular to p
        # start with a random vector
    n1 = np.random.rand(3)-0.5
        # subtract the component parallel to p to create a perpendicular vector 
    n1 -= np.dot(n1,p)*p
        # normalise to make it a unit vector
    n1 = n1/np.sqrt(np.dot(n1,n1))
    # make a second vector n2 that is perpendicular to both p and n1
    n2 = np.cross(p,n1)
    # n1 and n2 now span a plane n that is normal to p
    # find the point c where the plane n intersects the vector p for the given sin(alpha) sin_al
    cos_al = np.cos(np.arcsin(sin_al)) 
    # Scale the radius of the circle to sin(alpha)
    # produces a list of vectors that ascribe the circle
    ns = []
    for th in np.linspace(0, 2*np.pi, 1000):
        ns.append(sin_al*p+cos_al*(np.cos(th)*n1+np.sin(th)*n2))
    ns = np.array(ns)
    # check if they are unit vectors
    #print(np.sqrt(ns[:,0]**2+ns[:,1]**2+ ns[:,2]**2)) 
    # convert ns to azumuth/plunge (strike/dip) for plotting on stereonet
    strike, dip = mplstereonet.vector2pole(ns[:,0], ns[:,1], ns[:,2])
    return strike,dip
