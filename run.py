import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), 'fractoolbox'))


import fractoolbox as ftb

strike = ftb.dip2strike(45)
print(strike)