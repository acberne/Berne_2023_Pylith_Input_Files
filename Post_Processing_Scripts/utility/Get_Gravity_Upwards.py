

# *********************************************************************
# FUNCTION TO CONVERT SPHERICAL HARMONIC COEFFICIENTS INTO FROM A DENSITY INTERFACE INTO A GRAVITATIONAL ACCELERATION # #EVALUATED AT r, THETA, and PHI POINTS ABOVE INTERFACE at radius R up to SH degree lmax
#
# *********************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd
import pyshtools as pysh
from pyshtools import constants
from pyshtools import gravmag
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
import math




def main(r,thet,phi,lmax,clm,delrho,R):

    G=6.67e-11
    for l in range(1,lmax):
        for m in range(l+1):
            if m==0:
                norm=clm[0,l,m]
                vals[m,l,:]=((((norm)*delrho*(4*np.pi*G)*(pysh.expand.spharm_lm(l,m,theta*180/np.pi,phi*180/np.pi-((((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/l)))))*(R**(l+2))*(-(l+1))*(r**(-l-2))))*(1/(2*l+1))
            ##Next Step-- Sum up vals for all harmonics to that one point, Normalize, derivative wrong
                 ##-m
            else:
                norm=np.linalg.norm(clm[:,l,m])
                vals[m,l,:]=((((norm)*delrho*(4*np.pi*G)*(pysh.expand.spharm_lm(l,m,theta*180/np.pi,phi*180/np.pi-((((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/m)))))*(R**(l+2))*(-(l+1))*(r**(-l-2))))*(1/(2*l+1))		    


    return vals




