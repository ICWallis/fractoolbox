import numpy as np
import matplotlib.pyplot as plt
import functions as fun

# ============================================
# Model of thermal stress at the borehole wall
# ============================================
'''
Kirsh (1898) equations are used to calcuate the stresses resolved onto
the wall of a borehole where the axis is parallel to the vertical component
of the stress tensor. 
Equations as presented in Jager et al. (2007) and Zoback (2010).
'''

# Thermal Stress Magnatude
# ------------------------
'''
Calculate the magnatude of thermally induced stess [MPa] for a  
range of reservior tempratures and assuming steady state conditions.
The convention of - as tensile and + as compressive is used here.
The thermal stress is subtracted from the hoop stress curves 
and generate a shift in the negitive direction. 
'''
Twell = 50 + 273.14 # drilling gluid temp in Kelvin
therex = 1.e-5 # coefficent of thermal expansion
nu = 0.25 # Possions ratio
K = 1.e10 # bulk modulus

# thermal stress for reservior temps 50-300
rtemps = [50,100,150,200,250,300]
sigma_Dt_lst = []   
for n in rtemps:
    Tres = n + 273.14 # convert to Kelvin
    x = fun.fsigma_Dt(therex, K, nu, Tres, Twell)
    sigma_Dt_lst.append(x)

print(sigma_Dt_lst)

sigma_Dt = np.asarray(sigma_Dt_lst)
print(sigma_Dt)

# Geomechanical model paramaters
# ------------------------------
'''
Stress state and fluid pressures for a single depth 
'''

# strike slip faulting case
# stress/pressure in MPa
'''
SHmax = 90
Sv = 88.2
Shmin = 51.5 
Pp = 31.5  # why is this Pp so high?
Pmud = Pp 
'''

# normal faulting case
# stress/pressure in MPa
SHmax = 26
Sv = 88.232 
Shmin = 20 
Pp = 13.
Pmud = 14.75
deltaP = Pmud - Pp

# fixed paramaters 
nu = 0.1
R = 1 # wellbore radius
r = 1.0 # depth of investigation
# where R = r we are at the borehole wall
n = 200


plt.figure('Sigma rr',figsize=(10,5))

# angles around the borehole wall
theta = fun.ftheta(n)

# for f, b in zip(foo, bar):

for thermstress, rtemp in zip(sigma_Dt_lst,rtemps):
    tt = fun.effhoopstress(SHmax, Shmin, Pp, Pmud, thermstress, R, r, theta)
    plt.plot(theta*180/np.pi,tt,label=rtemp)


plt.xlabel(r'Angle around the borehole wall [$\theta$ degrees]')
plt.ylabel('Effective stress [MPa]')
plt.legend()
plt.show()

asdfag
# Radial stress
#rr = fun.sigma_rr(SHmax, Shmin, Pp, Pmud, R, r, theta)
#plt.plot(theta*180/np.pi,rr,label=r'$\sigma_{rr}$')

# Stress along the axis of the borehole
#zz = fun.effwbaxisstress(SHmax,Sv,Shmin,nu,Pp,R,r,theta,sigma_Dt) 
#plt.plot(theta*180/np.pi,zz,label=r'$\sigma_{zz}$')
# the mean of this should = the mean stress?






