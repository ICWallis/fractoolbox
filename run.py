import sys
import os

# sys.path.append(os.path.join(os.path.dirname(__file__), 'fractoolbox'))


from fractoolbox import dip2strike, thermal_stress

print(dip2strike(30))

print(thermal_stress(34, 2, 33, 22, 60))