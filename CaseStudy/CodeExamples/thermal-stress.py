# ============================================
# Model of thermal stress at the borehole wall
# ============================================
'''
Kirsh (1898) equations are used to calculate the stresses resolved onto the wall of 
a borehole. In this case, the borehole is vertical (i.e., the borehole axis is 
parallel to the vertical component of the stress tensor so the plots reflect the 
various influence of the two horizontal stresses as we move around the borehole wall. 
Equations as presented in Jager et al. (2007) and Zoback (2010).
'''
import numpy as np
import matplotlib.pyplot as plt
import functions as fun

# Thermal Stress Magnatude
# ------------------------
'''
Calculate the magnitude of thermally induced stress [MPa] for a  
range of reservoir temperatures and assuming steady state conditions.
The convention of - as tensile and + as compressive is used here.
The thermal stress is subtracted from the hoop stress curves 
and generate a shift in the negative direction. 
'''
Twell_degC = 50
Twell = Twell_degC + 273.14 # convert to kelvin
therex = 1.e-5 # coefficient of thermal expansion
nu = 0.25 # Possions ratio
K = 1.e10 # bulk modulus

# thermal stress for reservior temps 50-300
rtemps = [50,100,150,200,250,300]
thermal_stress_lst = []   
for n in rtemps:
    Tres = n + 273.14 # convert to Kelvin
    x = fun.thermal_stress(therex, K, nu, Tres, Twell)
    thermal_stress_lst.append(x)

# Model paramaters
# ----------------

# strike slip faulting case
# stress/pressure in MPa
SHmax_ss = 89
Sv_ss = 88
Shmin_ss = 40 

# normal faulting case
# stress/pressure in MPa
SHmax_nf = 26
Sv_nf = 88 
Shmin_nf = 20 

# fluid pressures
Pp = 13.
Pmud = 14.75

# fixed paramaters
R = 1 # wellbore radius
r = 1.0 # depth of investigation
# R = r at the borehole wall

# angles around the borehole wall
n = 200
theta = fun.theta(n)

# Calculate and plot
# ------------------
f, (ax1,ax2) = plt.subplots(1, 2, figsize = (15, 5), sharey=True)

# hoop stress - normal faulting case
for thermstress, rtemp in zip(thermal_stress_lst,rtemps):
    tt = fun.effhoopstress(SHmax_nf, Shmin_nf, Pp, Pmud, thermstress, R, r, theta)
    ax1.plot(theta*180/np.pi,tt,label=rtemp)

# hoop stress - strike slip case
for thermstress, rtemp in zip(thermal_stress_lst,rtemps):
    tt = fun.effhoopstress(SHmax_ss, Shmin_ss, Pp, Pmud, thermstress, R, r, theta)
    ax2.plot(theta*180/np.pi,tt,label=rtemp)

# tensile rock strength
for ax in [ax1,ax2]: 
    ax.hlines(1, 0, 360,
          colors='k', 
          linestyles='--', 
          linewidth=2,
          alpha=0.2)

# compressive rock strength
UCS = [(20,30),(70,100)]
for ax in [ax1,ax2]:
    for top, bottom in UCS:
        ax.axhspan(top, bottom, color='k', alpha=.2)
'''
From Wyering et al (2014)
-   shallower samples with low-temperature alteration UCS = 27.7 ± 10.3 MPa
-   deeper samples with high-temperature alteration 84.8 ± 30.6 MPa
Range used here is 20-30 for shallow and 70-100 for deep
'''

for ax in [ax1,ax2]:
    ax.set_xlabel(r'Angle on borehole wall [deg measured from the SHmax azumuth]')
    ax.set_xlim(0,360)
    ax.set_xticks([0,90,180,270,360])
    ax.legend(loc='upper right')

ax1.set_ylim(-60,120)
ax1.set_ylabel('Effective stress [MPa]')

ax1.set_title('Normal faluting (small difference in horizontal stress)' )
ax2.set_title('Strike-slip faluting (large difference in horizontal stress)' )

plt.suptitle('Same depth (Sv) with a fixed ' + str(Twell_degC) + 'degC mud and varying reservior temp')
plt.show()