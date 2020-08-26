# this method has been coded from PP 156-157 of Zoback 2010

datapath = r'/Users/irene/Dropbox/PhD/WellData/ProcessedData/'
figexportpath = r'/Users/irene/Dropbox/PhD/Manuscripts/NZGW20 Feedzones/BackgroundFigures/'
# Python version = 3.6
# conda env = tank

from matplotlib import pyplot as plt
import numpy as np
import math
import pandas as pd
import mplstereonet

pd.set_option('display.width', 1000)
pd.set_option('max_colwidth', 100) # made this longer because my tuple is so long that it will not display

# my library
#import stress as sts
#import wsgtools as wsg
import functions as fun

# Jobs to do
# ----------
# does the coloumb falure criterion get normalised to effective vertical stress somehow? It seems to me to be just a gradient.
# go back to the normalisation issue, because we are not generating ovelapping mohr plots again

# -----------------------
# IMPORT AND MASSAGE DATA
# -----------------------

# test fractures
# -------------------------------------------------------------------------------
# 
dffracs = pd.read_csv(r'demo-presentation/3Dmohr-plot-fractures.csv')

# append stress field orentation and make a tuple to pass into the function
for df in [dffracs]:
    # alpha = -60, beta = -60, gamma = 0 gives the best look for ignimbrite only
    df['alpha'] = -60 
    df['beta'] = -90 # -90 is the origional and -60 is the tilted
    df['gamma'] = 0 
    df['Sv_eff'] = df.Sv - df.Pp
    df['SHmax_eff'] = df.SHmax - df.Pp
    df['Shmin_eff'] = df.Shmin - df.Pp 
    # so I am using effective vertical stress rather than vertical stress without changing the functions
    # Sv in the fucntion list below is just used to normalise... perhaps I should change the name?
    # compile my data into a tuple that can be used in the function using the following sequence (S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip)
    df['fracture'] = df[['Sv','SHmax','Shmin','Pp','Sv_eff','alpha','beta','gamma','strike','dip']].apply(tuple, axis=1)

# ----------
# THE METHOD
# ----------
# Do not change anything in this section
# At some point this should be moved to the stress function library, 
# but I would need to make the object definitons unique to do so

def FractureSnTau(S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip):
    '''Calculate the shear (tau) and normal (Sn) stress on a fracture

    NOTE I need to check if Sv should really be Sv effective 
    (because I've normalised by effective stress elsewhere)

    Args:
        S1:
        S2:
        S3:
        Pp:
        Sv:
        alpha:
        beta:
        gamma:
        strike:
        dip:

        Reccomendation: this can be efficently done using a tuple 
    
    Returns: 
        Sn: Stress normal to the fracture plane [MPa]
        tau: Shear stress on the fracture plane [MPa]

    '''
    # create effective stress array
    Ss = np.array([                 
                    [S1-Pp,0,0],
                    [0,S2-Pp,0],
                    [0,0,S3-Pp]
                    ])
    #print('Ss: the effective stress tensor =','\n',Ss,'\n')

    # use the three Euiler angles to generate an array that is used 
    # to transform the effective stress array into geographic coordinates
    Rs = fun.Rs(alpha,beta,gamma)  
    #print('Rs: the stress coordinate system based on' + 
    #   'the inputted Euler angles =','\n',Rs,'\n')

    # rotate the stress tensor into geographic cooridnates
    Sg = Rs.T@Ss@Rs                 
    #print('Sg: the effective stress tensor now rotated into' +
    #   'geographic coordinates =','\n',Sg,'\n')

    # use the fracture strike an dip to generate an array that is used
    # to transform the stress field into fracture cooridinates
    Rf = fun.Rf(strike,dip)        
    #print('Rf: the matrix that rotates the stress tensor from' + '
    #   'geographic coordinates into the fracture plane coordinates =','\n',Rf,'\n')

    # transform the stress field into the fracture coordinate system
    Sf = Rf@Sg@Rf.T                 
    #print('Sf: the effective stress tensor now rotated into the' +'
    #   fracture plane coordinate system =','\n',Sf,'\n')

    # take stress normal to the fault plane and normalise it to 
    # vertical stress so fractures from different depths can be plotted together
    Sn = Sf[2,2]/Sv                 
    #print('Sn: the effective stress magnitude normal to the' + '
    #   'fault plane (Sf bottom right) normalised to Sv =','\n',Sn,'\n')

    # calcuate the rake of the fracture assuming only dip-slip
    rake = fun.rake(Sf)            
    #print('the rake of the slip vector =',rake,'\n')

    # use the rake to generate an array to transform the stress field into the rake
    Rt = fun.Rt(rake)
    #print('Array Rt that is used to transform the effective stress from the' + 
    #   'fault co-ordinate system into the rake coordinate system so we can' + 
    #   'grab the shear stress magnitude along the slip vector =','\n',Rt,'\n')

    # transform the stress field into the direction of the rake
    Sr = Rt@Sf@Rt.T
    #print('Effective stress tensor along the slip vector (Sr) where the' + 
    #   'bottom left (as an absolute number) is the shear stress magnitude (tau) =','\n',Sr,'\n')

    # take the absolue number of shear stress in the direction of the rake 
    # and normalise to vertical stress for plotting
    tau = abs(Sr[2,0])/Sv
    #print('Shear stress on the plane (tau, which is bottom left of Sr array) =',tau,'\n')

    return (Sn,tau)

# ---------------------------------
# CALCULATE SHEAR AND NORMAL STRESS
# ---------------------------------
# Calls the method above to calculate the shear and normal stress for each fracture
# method is set up to loop over multiple df when required

for df in [dffracs]:
    dalist = []
    for  S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip in df['fracture']:
        Sn, tau = FractureSnTau(S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip)
        dalist.append([Sn, tau])
    x = pd.Series(dalist)
    df['Sn_tau'] = x.values


# ---------------------------------------------
# CALCULATE THE RATIO OF SHEAR TO NORMAL STRESS
# ---------------------------------------------
# method is set up to loop over multiple df

for df in [dffracs]:
    dalist = []
    for Sn, tau, in df['Sn_tau']:
        ratio = tau/Sn
        dalist.append(ratio)
    x = pd.Series(dalist)
    df['ratio'] = x.values

print(dffracs)
 

# GENERATE MOHR PLOT OUTLINE
# --------------------------
# perhaps the problem is because i am taking absolute magnatudes rather than the magnatudes in geographic cooridnates?
# because in geographic coordinates, I lose a bit in the rotation?

# calling the shallowest (smallest) effective stress values
# ---------------------------------------------------------

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

# Rearrange the Mohr values for plotting
# --------------------------------------
# shallow stress values
pltpairs = {                        
    "normStauS":
    [(normS[0],tauS[0]),            # sigma1, sigma3 circle
    (normS[1],tauS[1]),             # sigma1, sigma2 circle
    (normS[2],tauS[2])]             # sigma2, sigma3 circle
}

# -----------------------
# COLOMB FALURE CRITERION
# -----------------------
# Generate data to make the criterion lines on the plots 
# Following Barton et al 1995

CFCSn = [0,100]
CFCtau_pnt5 = []
for n in CFCSn:
    x = 0.5 * n                     # mu of 0.5 - this is what we used in the XLOT paper because it matched the measurements
    CFCtau_pnt5.append(x)

CFCtau_pnt6 = []
for n in CFCSn:
    x = 0.6 * n                     # mu of 0.6
    CFCtau_pnt6.append(x)

CFCtau_1 = []
for n in CFCSn:
    x = 1 * n                       # mu of 1
    CFCtau_1.append(x)


# -------------
# DRAW THE PLOT
# -------------

f = plt.figure(figsize=(8,4))
ax1 = f.add_subplot(111)

# plot the coloumb falure criterion
ax1.plot(CFCSn,CFCtau_1,c='k',linewidth=0.5)
ax1.plot(CFCSn,CFCtau_pnt6,c='k',linewidth=0.5)
ax1.plot(CFCSn,CFCtau_pnt5,c='k',linewidth=0.5) 
ax1.text(0.35, 0.4, r'$\mu = 1$', fontsize=10, bbox={'facecolor': 'white','edgecolor':'none', 'pad': 2})
ax1.text(0.57, 0.4, r'$\mu = 0.6$', fontsize=10, bbox={'facecolor': 'white','edgecolor':'none', 'pad': 2})
ax1.text(0.77, 0.4, r'$\mu = 0.5$', fontsize=10, bbox={'facecolor': 'white','edgecolor':'none', 'pad': 2})

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