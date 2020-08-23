'''
============================
fractoolbox Function Library
============================



'''

# Convert Fracture Data Format
# ============================
'''

'''

def dip2strike(dipazim,dip):
    '''Convert dip-dipazimuth data to strike 
    using the right-hand rule convention'''
    if dipazim < 90:
        strike=(dipazim-90)+360
    else:
        strike=dipazim-90
    return strike

# Geometric Sample Bias: Isogenic Contours
# ========================================
'''

'''