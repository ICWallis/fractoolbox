# Coded from PP 156-157 of Zoback 2010
# Python version = 3.6
# conda env = tank
# Normalising to effective stress throughout

from matplotlib import pyplot as plt
import numpy as np
import math
import pandas as pd
import mplstereonet

import functions as fun

pd.set_option('display.width', 1000)
# extend to tuple of params will display
pd.set_option('max_colwidth', 100) 

# Import data and setup dataframe
# --------------------------------
dffracs = pd.read_csv(r'demo-presentation/3Dmohr-plot-fractures.csv')

# add three fractures to the dataframe, one each that is perpendicular to the principal stresses
# the stress and pressure values come from the synthetic fracture at 2000 mVD
#print(list(dffracs.columns))
# order = ['ID', 'mVD', 'dip', 'dipAZ', 'strike', 'Pp', 'Sv', 'Shmin', 'SHmax']
# note that dip az position is empty because it's not used in this method
extras = [
        pd.Series([5, 2000, 0.0001, '', 30, 16.41, 39.20, 25.12, 32.16], index=dffracs.columns ), 
        pd.Series([6, 2000, 90, '', 30, 16.41, 39.20, 25.12, 32.16], index=dffracs.columns),
        pd.Series([7, 2000, 90, '', 120, 16.41, 39.20, 25.12, 32.16], index=dffracs.columns) 
        ]
dffracs = dffracs.append(extras, ignore_index=True)

# Synthetic data for play example
# -------------------------------
'''
fracture = {
    'ID': [],
    'mVD': [],
    'dip': [],
    'strike': [],
    'Pp': [],
    'Sv': [],
    'Shmin': [],
    'SHmax': [],
}

dftest = pd.DataFrame(data=fracture)
dftest
'''
# append stress field orentation and make a tuple for Sn Tau calculation
# ----------------------------------------------------------------------
# alpha: az of SHmax
# beta: tils of Sv, -90 is the origional and -60 is the 30 deg tilted for the case study
for df in [dffracs]:
    df['alpha'] = -60 
    df['beta'] = -90 
    df['gamma'] = 0 
    df['Sv_eff'] = df.Sv - df.Pp
    df['SHmax_eff'] = df.SHmax - df.Pp
    df['Shmin_eff'] = df.Shmin - df.Pp 
    # make tuple for fracture_sn_tau(S1,S2,S3,Pp,Norm,alpha,beta,gamma,strike,dip)
    df['fracture'] = df[['Sv','SHmax','Shmin','Pp','Sv_eff','alpha','beta','gamma','strike','dip']].apply(tuple, axis=1)

# Calculate shear and normal stress magnatude
# -------------------------------------------
for df in [dffracs]:
    dalist = []
    for  S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip in df['fracture']:
        Sn, tau = fun.fracture_sn_tau(S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip)
        dalist.append([Sn, tau])
    x = pd.Series(dalist)
    df['Sn_tau'] = x.values

# Calculate the ratio of shear to normal stres
# --------------------------------------------
for df in [dffracs]:
    dalist = []
    for Sn, tau, in df['Sn_tau']:
        ratio = tau/Sn
        dalist.append(ratio)
    x = pd.Series(dalist)
    df['ratio'] = x.values

# Generate Mohr plot parts
# --------------------------

# Min Sv value (use effective stress)
Sigma1 = dffracs['Sv_eff'].min()
print(Sigma1)

# Min SHmax value
Sigma2 = dffracs['SHmax_eff'].min()
print(Sigma2)

# Min Shmin value
Sigma3 = dffracs['Shmin_eff'].min()
print(Sigma3)

# make the components for the 3D mohr plot
tauS, normS, meanS = fun.mohr3d(Sigma1,Sigma2,Sigma3,Sigma1) 

# Rearrange for plotting
pltpairs = {                        
    "normStauS":
    [(normS[0],tauS[0]),  # sigma1, sigma3 circle
    (normS[1],tauS[1]),   # sigma1, sigma2 circle
    (normS[2],tauS[2])]   # sigma2, sigma3 circle
}

# Make coloumb falure criterion
# -----------------------------
# Generate data to make the criterion lines on the plots 
# following the format used in Barton et al 1995

CFCSn = [0,100]

# mu of 0.5
CFCtau_pnt5 = []
for n in CFCSn:
    x = 0.5 * n
    CFCtau_pnt5.append(x)

# mu of 0.6
CFCtau_pnt6 = []
for n in CFCSn:
    x = 0.6 * n
    CFCtau_pnt6.append(x)

# mu of 1
CFCtau_1 = []
for n in CFCSn:
    x = 1 * n
    CFCtau_1.append(x)

# Draw the plot
# -------------
f = plt.figure(figsize=(8,4))
ax1 = f.add_subplot(111)

# plot the coloumb falure criterion
ax1.plot(CFCSn,CFCtau_1,c='k',linewidth=0.5)
ax1.plot(CFCSn,CFCtau_pnt6,c='k',linewidth=0.5)
ax1.plot(CFCSn,CFCtau_pnt5,c='k',linewidth=0.5) 
ax1.text(0.35, 0.4, r'$\mu = 1$', fontsize=10, 
    bbox={'facecolor': 'white','edgecolor':'none', 'pad': 2})
ax1.text(0.57, 0.4, r'$\mu = 0.6$', fontsize=10, 
    bbox={'facecolor': 'white','edgecolor':'none', 'pad': 2})
ax1.text(0.77, 0.4, r'$\mu = 0.5$', fontsize=10, 
    bbox={'facecolor': 'white','edgecolor':'none', 'pad': 2})

# plot the mean stress tick marks
ax1.scatter(meanS,[0,0,0],color='k',marker='|')

# plot the Mohr plot arcs
for normS, tauS in pltpairs['normStauS']:
    ax1.plot(normS,tauS,'k',linewidth=0.5)

# Plot the fractures
for Sn, tau, in dffracs['Sn_tau']:
    ax1.scatter(Sn,tau,s=6,c='k',marker='o')

# Final formatting and output
ax1.set_ylim(0,0.45)
ax1.set_xlim(0,1.2)
ax1.set_aspect('equal')
ax1.grid(linestyle='--', linewidth=0.2)

ax1.set_xlabel(r'$(\sigma_n - Pp)/(Sv-Pp)$')
ax1.set_ylabel(r'$ \tau_s/(Sv-Pp)$')

#f.subplots_adjust(top=10)
#plt.suptitle('title',fontsize=12)

plt.show()