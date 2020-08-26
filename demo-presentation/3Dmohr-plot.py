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
    Rs = sts.fRs(alpha,beta,gamma)  
    #print('Rs: the stress coordinate system based on the inputted Euler angles =','\n',Rs,'\n')

    # rotate the stress tensor into geographic cooridnates
    Sg = Rs.T@Ss@Rs                 
    #print('Sg: the effective stress tensor now rotated into geographic coordinates =','\n',Sg,'\n')

    # use the fracture strike an dip to generate an array that is used
    # to transform the stress field into fracture cooridinates
    Rf = sts.fRf(strike,dip)        
    #print('Rf: the matrix that rotates the stress tensor from geographic coordinates into the fracture plane coordinates =','\n',Rf,'\n')

    # transform the stress field into the fracture coordinate system
    Sf = Rf@Sg@Rf.T                 
    #print('Sf: the effective stress tensor now rotated into the fracture plane coordinate system =','\n',Sf,'\n')

    Sn = Sf[2,2]/Sv                 # take stress normal to the fault plane and normalise it to vertical stress for plotting
    #print('Sn: the effective stress magnitude normal to the fault plane (Sf bottom right) normalised to Sv =','\n',Sn,'\n')

    rake = sts.frake(Sf)            # calcuate the rake of the fracture
    #print('the rake of the slip vector =',rake,'\n')

    Rt = sts.fRt(rake)               # use the rake to generate an array to transform the stress field into the rake
    #print('Array Rt that is used to transform the effective stress from the fault co-ordinate system into the rake coordinate system so we can grab the shear stress magnitude along the slip vector =','\n',Rt,'\n')

    Sr = Rt@Sf@Rt.T                 # transform the stress field into the direction of the rake
    #print('Effective stress tensor along the slip vector (Sr) where the bottom left (as an absolute number) is the shear stress magnitude (tau) =','\n',Sr,'\n')

    tau = abs(Sr[2,0])/Sv           # take the absolue number of shear stress in the direction of the rake for plotting
    #print('Shear stress on the plane (tau, which is bottom left of Sr array) =',tau,'\n')

    return (Sn,tau)

# ---------------------------------
# CALCULATE SHEAR AND NORMAL STRESS
# ---------------------------------
# Calls the method above to calculate the shear and normal stress for each fracture
# Use one if the loops below and comment out the other

# Loop for when there is one dataframe
# ------------------------------------
# eg the synthetic fracture dataframe
'''
dalist = []
for  S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip in df['fracture']:
    Sn, tau = FractureSnTau(S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip)
    dalist.append([Sn, tau])
x = pd.Series(dalist)
df['Sn_tau'] = x.values
'''
# Loop for when there are multiple dataframes
# -------------------------------------------
# eg the NM10 dataset filtered by lithology

for df in [dfconductivehalo,dfconductivenonhalo,dffaults,dfresistive,dfiso_0pnt9,dfiso_0pnt4,dfiso_0pnt3,dfiso_0pnt2,dfiso_0pnt1]:
    dalist = []
    for  S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip in df['fracture']:
        Sn, tau = FractureSnTau(S1,S2,S3,Pp,Sv,alpha,beta,gamma,strike,dip)
        dalist.append([Sn, tau])
    x = pd.Series(dalist)
    df['Sn_tau'] = x.values

# ---------------------------------------------
# CALCULATE THE RATIO OF SHEAR TO NORMAL STRESS
# ---------------------------------------------
# Use one if the loops below and comment out the other

# Loop for when there is one dataframe
# ------------------------------------
# eg the synthetic fracture dataframe
'''
dalist = []
for Sn, tau, in df['Sn_tau']:
    ratio = tau/Sn
    dalist.append(ratio)
x = pd.Series(dalist)
df['ratio'] = x.values
'''
# Loop for when there is one dataframe
# ------------------------------------
# eg the synthetic fracture dataframe

for df in [dfconductivehalo,dfconductivenonhalo,dffaults,dfresistive,dfiso_0pnt9,dfiso_0pnt4,dfiso_0pnt3,dfiso_0pnt2,dfiso_0pnt1]:
    dalist = []
    for Sn, tau, in df['Sn_tau']:
        ratio = tau/Sn
        dalist.append(ratio)
    x = pd.Series(dalist)
    df['ratio'] = x.values

# Investigation or filter the ratios
# ----------------------------------
#print(dfsed['ratio'])

# Make a histogram of the ratio
#n, bins, patches = plt.hist(dfand['ratio'],50,density=True,facecolor = 'g', alpha=0.75) 

# Make a new dataframe containing those fractures above a set value
#dfandPnt5 = dfand.loc[dfand.ratio > 0.5, :] 

# GENERATE MOHR PLOT OUTLINE
# --------------------------
# perhaps the problem is because i am taking absolute magnatudes rather than the magnatudes in geographic cooridnates?
# because in geographic coordinates, I lose a bit in the rotation?

# calling the shallowest (smallest) stress values
# -----------------------------------------------
# for ease today, I read these values straight from the dfs in the data processing ipynb
# these convert the stresses to effective stresses 

### NOTE ### CHECK ### NOTE ### UPDATE ### NOTE ###
### NOTE ### CHECK ### NOTE ### UPDATE ### NOTE ###
### NOTE ### CHECK ### NOTE ### UPDATE ### NOTE ###

PoreP = 16.06
print(Pp)

# Min Sv value
Sigma1 = 38.90 - PoreP 
print(Sigma1)

# Min SHmax value
Sigma2 = 31.14  - PoreP
print(Sigma2)

# Min Shmin value
Sigma3 = 23.38 - PoreP
print(Sigma3)

# call the function that makes all of the components for the mohr plot
# remeber that the last object passed into this is the normalisation value
tauS, normS, meanS = sts.fMohr3DSvnorm(Sigma1,Sigma2,Sigma3,Sigma1) 

# It is not yet clear if I need to make this rotation
# MORE WORK REQUIRED
'''
# making the rotation and then drawing the values... do we need to rotate?? no??
S = sts.fSsEf_GeoCoords(Sigma1,Sigma2,Sigma3,PoreP,30,0,90) # make the roation into geographic co-oridates for this one case
print(S)

# array positions for a strike slip case
S1 = S[0,0] # S1
S2 = S[2,2] # S2
S3 = S[1,1] # S3
N = S[2,2]  # bottom right is the vertical stress in all cases

tauS, normS, meanS = sts.fMohr3DSvnorm(S1,S2,S3,N) # call the function that makes all of the components for the mohr plot
'''

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

# Construct the 3D Mohr circles
# -----------------------------
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

# Plot the isogenic contours
# --------------------------

#for Sn, tau, in dfiso_0pnt1['Sn_tau']:
#    ax1.scatter(Sn,tau,s=3,c='#eee9f8')

#for Sn, tau, in dfiso_0pnt2['Sn_tau']:
#    ax1.scatter(Sn,tau,s=3,c='#cdbeec')

#for Sn, tau, in dfiso_0pnt3['Sn_tau']:
#    ax1.scatter(Sn,tau,s=3,c='#ac93e0')

#for Sn, tau, in dfiso_0pnt4['Sn_tau']:
#    ax1.scatter(Sn,tau,s=3,c='#8b67d3')

#for Sn, tau, in dfiso_0pnt9['Sn_tau']:
#    ax1.scatter(Sn,tau,s=3,c='#3e1b87')

# Plot the fractures
# ------------------
# Filter the dataframes by lithology 
# Refer to 3DMohrPlot-ProcessLogFractureDataset-Normal.ipynb Step 3 to generate and test these

# Andesite (all occurance types in the stratigraphic unit)

dfconductivenonhalo_andesite = dfconductivenonhalo[
    (dfconductivenonhalo.Lith_PetroFacies == 'Lower Andesite') |
    (dfconductivenonhalo.Lith_PetroFacies == 'Upper Andesite') |
    (dfconductivenonhalo.Lith_PetroFacies == 'Other (Intrusive?)')
].copy()

dfconductivehalo_andesite = dfconductivehalo[
    (dfconductivehalo.Lith_PetroFacies == 'Lower Andesite') |
    (dfconductivehalo.Lith_PetroFacies == 'Upper Andesite') |
    (dfconductivehalo.Lith_PetroFacies == 'Other (Intrusive?)')
].copy()

dffaults_andesite = dffaults[
    (dffaults.Lith_PetroFacies == 'Lower Andesite')
].copy()

dfresistive_andesite = dfresistive[
    (dfresistive.Lith_PetroFacies == 'Lower Andesite') |
    (dfresistive.Lith_PetroFacies == 'Upper Andesite') |
    (dfresistive.Lith_PetroFacies == 'Other (Intrusive?)')
].copy()

# all clastic units

dfconductivenonhalo_clastic = dfconductivenonhalo[
    (dfconductivenonhalo.Lith_PetroFacies == 'Mixed Volcaniclastics')
].copy()

dfconductivehalo_clastic = dfconductivehalo[
    (dfconductivehalo.Lith_PetroFacies == 'Mixed Volcaniclastics')
].copy()

dffaults_clastic = dffaults[
    (dffaults.Lith_PetroFacies == 'Mixed Volcaniclastics')
].copy()

dfresistive_clastic = dfresistive[
    (dfresistive.Lith_PetroFacies == 'Mixed Volcaniclastics')
].copy()

# ignimbrite including the ash layers

dfconductivenonhalo_ignimbrite = dfconductivenonhalo[
    (dfconductivenonhalo.Lith_PetroFacies == 'Pyroclastic')
].copy()

dfconductivehalo_ignimbrite = dfconductivehalo[
    (dfconductivehalo.Lith_PetroFacies == 'Pyroclastic')
].copy()

dffaults_ignimbrite = dffaults[
    (dffaults.Lith_PetroFacies == 'Pyroclastic') |
    (dffaults.Lith_PetroFacies == 'Ash')
].copy()

dfresistive_ignimbrite = dfresistive[
    (dfresistive.Lith_PetroFacies == 'Pyroclastic')
].copy()

# Plot the andesite filtered data
# -------------------------------
for Sn, tau, in dfconductivenonhalo_andesite['Sn_tau']:
    ax1.scatter(Sn,tau,s=6,c='#7f7f7f',marker='o')  # use #b2d1ce if plotting with the other colours

#for Sn, tau, in dfconductivehalo_andesite['Sn_tau']:
#    ax1.scatter(Sn,tau,s=10,c='#1a756e',marker='^')

#for Sn, tau, in dffaults_andesite['Sn_tau']:
#    ax1.scatter(Sn,tau,s=10,c='#00332f',marker='D')

#for Sn, tau, in dfresistive_andesite['Sn_tau']:
#    ax1.scatter(Sn,tau,s=8,c='#1a756e',marker='s')

# Plot the ignimbrite filtered data
# ---------------------------------
#for Sn, tau, in dfconductivenonhalo_ignimbrite['Sn_tau']:
#    ax1.scatter(Sn,tau,s=6,c='#7f7f7f',marker='o')  # use #d1b99d if plotting with the other colours

#for Sn, tau, in dfconductivehalo_ignimbrite['Sn_tau']:
#    ax1.scatter(Sn,tau,s=10,c='#976222',marker='^')

#for Sn, tau, in dffaults_ignimbrite['Sn_tau']:
#    ax1.scatter(Sn,tau,s=10,c='#462805',marker='D')

#for Sn, tau, in dfresistive_ignimbrite['Sn_tau']:
#    ax1.scatter(Sn,tau,s=8,c='#976222',marker='s')

# Plot the clastic filtered data (same schema as above)
# -----------------------------------------------------
#for Sn, tau, in dfconductivenonhalo_clastic['Sn_tau']:
#    ax1.scatter(Sn,tau,s=6,c='#7f7f7f',marker='o')  # use #d1b99d if plotting with the other colours

#for Sn, tau, in dfconductivehalo_clastic['Sn_tau']:
#    ax1.scatter(Sn,tau,s=10,c='#976222',marker='^')

#for Sn, tau, in dffaults_clastic['Sn_tau']:
#    ax1.scatter(Sn,tau,s=10,c='#462805',marker='D')

#for Sn, tau, in dfresistive_clastic['Sn_tau']:
#    ax1.scatter(Sn,tau,s=8,c='#976222',marker='s')

# all conductive non-halo
# -----------------------
#for Sn, tau, in dfconductivenonhalo['Sn_tau']:
#    ax1.scatter(Sn,tau,s=4,c='#7f7f7f',marker='o')


# Final formatting and output
# ---------------------------
ax1.set_ylim(0,0.45)
ax1.set_xlim(0,1.2)
ax1.set_aspect('equal')
ax1.grid(linestyle='--', linewidth=0.2)

ax1.set_xlabel(r'$(\sigma_n - Pp)/(Sv-Pp)$')
ax1.set_ylabel(r'$ \tau_s/(Sv-Pp)$')

#f.subplots_adjust(top=10)
#plt.suptitle('Mohr circle testing with a normal faulting case (alpha = -60, beta = -90, gamma = 0)',fontsize=12)

#plt.savefig(figexportpath + r'3DMohrPlots-NM10-plot-normal-Andesite.png',dpi=300,bbox_inches='tight')
#plt.savefig(figexportpath + r'3DMohrPlots-NM10-plot-normal-IgAndClastic.png',dpi=300,bbox_inches='tight')
#plt.savefig(figexportpath + r'3DMohrPlots-NM10-plot-normal-AllConductive.png',dpi=300,bbox_inches='tight')

#plt.savefig(figexportpath + r'3DMohrPlots-NM10-plot-normal-Ig-Tilted.png',dpi=300,bbox_inches='tight')
#plt.savefig(figexportpath + r'3DMohrPlots-NM10-plot-normal-Ig-Andersonian.png',dpi=300,bbox_inches='tight')
#plt.savefig(figexportpath + r'3DMohrPlots-NM10-plot-normal-Ande-Andersonian.png',dpi=300,bbox_inches='tight')
#plt.savefig(figexportpath + r'3DMohrPlots-NM10-plot-normal-Ande-Tilted.png',dpi=300,bbox_inches='tight')


# export df for use in stereonet plots
# one for each of the stress orentation cases
#dfconductivenonhalo.to_csv(datapath + r'NM10-Fracs-SLB-ConductiveFracture-NonHalo-Merged-Snf-alpha-60beta-90gamma0.csv',index = None, header=True)
#dfconductivenonhalo.to_csv(datapath + r'NM10-Fracs-SLB-ConductiveFracture-NonHalo-Merged-Snf-alpha-60beta-60gamma0.csv',index = None, header=True)

#plt.show()