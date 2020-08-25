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
Twell_degC = 50
Twell = Twell_degC + 273.14 # convert to kelvin
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

f, ax = plt.subplots(1, 1, figsize = (10, 5))

#plt.figure('Sigma rr',figsize=(10,5))

# angles around the borehole wall
theta = fun.ftheta(n)

for thermstress, rtemp in zip(sigma_Dt_lst,rtemps):
    tt = fun.effhoopstress(SHmax, Shmin, Pp, Pmud, thermstress, R, r, theta)
    ax.plot(theta*180/np.pi,tt,label=rtemp)

# tensile rock strength
ax.hlines(1,
          0,360,
          colors='k', 
          linestyles='--', 
          linewidth=2,
          alpha=0.2)

# compressive rock strength
'''
From Wyering et al (2014)
-   shallower samples with low-temperature alteration UCS = 27.7 ± 10.3 MPa
-   deeper samples with high-temperature alteration 84.8 ± 30.6 MPa
Range used here is 20-30 for shallow and 70-100 for deep
'''
UCS = [(20,30),(70,100)]

for top, bottom in UCS:
    ax.axhspan(top, bottom, color='k', alpha=.2)

ax.set_xlabel(r'Angle around the borehole wall [$\theta$ degrees]')
ax.set_ylabel('Effective stress [MPa]')
ax.set_xlim(0,360)
ax.set_ylim(-60,120)
ax.set_xticks([0,90,180,270,360])
ax.legend(loc='upper right')

ax.set_title('Stress resolved onto the borehole wall, with fixed ' + str(Twell_degC) + 'degC mud and varying reservior temp' )
plt.show()

# Radial stress
#rr = fun.sigma_rr(SHmax, Shmin, Pp, Pmud, R, r, theta)
#plt.plot(theta*180/np.pi,rr,label=r'$\sigma_{rr}$')

# Stress along the axis of the borehole
#zz = fun.effwbaxisstress(SHmax,Sv,Shmin,nu,Pp,R,r,theta,sigma_Dt) 
#plt.plot(theta*180/np.pi,zz,label=r'$\sigma_{zz}$')
# the mean of this should = the mean stress?






