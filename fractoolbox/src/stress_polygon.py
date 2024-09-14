'''In progress

'''
import numpy as np
import math
import pandas as pd
from scipy import integrate
from matplotlib import pyplot as plt


def minstress(S1,Pp,mu):
    '''Use the stress ratio and frictional faulting theory to estimate the minimum stress
    
    Args:
        Pp: Pore pressure MPa 
        S1: Maximum stress MPa
        
    Returns:
        S3: Minimum stress    
    
    '''
    S3 = ((S1-Pp)/(((mu**2 + 1.)**0.5 + mu)**2))+Pp
    return S3

def maxstress(S3,Pp,mu):
    '''Using the stress ratio and frictional faulting theory to estimate the maximum stress
    
    Args:
        Pp: Pore pressure MPa 
        S3: Minimum stress MPa
        
    Returns:
        S1: Maximum stress MPa
        
    '''
    S1 = ((S3-Pp)*(((mu**2 + 1.)**0.5 + mu)**2))+Pp
    return S1

def poly(Sv,Pp,mu,figname='StressPolygon'):
    '''Draws a stress polygon plot

    The stress polygon is the minimum and maximum horizontal stresses allowable based on the 
    Mohr-coloumb failure criterion. The edge of the polygon is where failure on an 
    optimally orientated fault or fracture will occur. 
    
    Args:
        Sv: Vertical stress [MPa]
        Pp: Pore pressure [MPa]
        mu: Coefficient of friction (0.1-0.6 with 0.5 a reasonable first estimate)

    Returns:
        A plot of stress polygons sometimes referred to as the Zoback-a-gram    
    '''
    minSh = minstress(Sv,Pp,mu)
    maxSh = maxstress(Sv,Pp,mu)
    minSH = minSh
    maxSH = maxSh
    ax = [minSh,minSh] # endpoints of the connecting lines
    ay = [minSh,Sv]
    bx = [minSh,Sv]
    by = [Sv,maxSH]
    cx = [Sv,maxSh]
    cy = [maxSH,maxSH]
    dx = [minSh,Sv]
    dy = [Sv,Sv]
    ex = [Sv,Sv]
    ey = [Sv,maxSH]
    fx = [minSh,maxSH]
    fy = [minSh,maxSH]
    f,ax1 = plt.subplots(1,1,figsize=(6,6))
    ax1.plot(ax,ay,color='k',alpha=0.5) # plots the connecting lines
    ax1.plot(bx,by,color='k',alpha=0.5)
    ax1.plot(cx,cy,color='k',alpha=0.5)
    ax1.plot(dx,dy,color='k',alpha=0.5)
    ax1.plot(ex,ey,color='k',alpha=0.5)
    ax1.plot(fx,fy,color='k',alpha=0.5)
    ax1.plot(minSh,minSH,'o',color='k')  # 1. Highly extensional $ S_hmin = S_Hmax << S_v $
    ax1.plot(Sv,Sv,'o',color='k')        # 2. Central point $ S_hmin = S_Hmax = S_v $
    ax1.plot(minSh,Sv,'o',color='k')     # 3. Transition between NF and SS $ S_hmin < S_Hmax = S_v $
    ax1.plot(Sv,maxSH,'o',color='k')     # 4. Transition between SS and RF $ S_hmin = S_v << S_Hmax$
    ax1.plot(maxSh,maxSH,'o',color='k')  # 5. Highy compresional $ S_hmin = S_Hmax >> S_v $
    ax1.text((maxSH-Sv)/4+Sv,(maxSH-Sv)/1.5+Sv, 'RF', fontsize=10)
    ax1.text((Sv-minSh)/2+minSh,(maxSH-Sv)/8+Sv, 'SS', fontsize=10)
    ax1.text((Sv-minSh)/4+minSh,(Sv-minSh)/1.5+minSh, 'NF', fontsize=10)
    plt.xlim(0,maxSh+20)
    plt.ylim(0,maxSh+20)
    plt.grid(linestyle='--')
    plt.xlabel('$S_{hmin}$ [MPa]')
    plt.ylabel('$S_{Hmax}$ [MPa]')
    plt.savefig(figname + '.png', dpi=300)
    plt.show()