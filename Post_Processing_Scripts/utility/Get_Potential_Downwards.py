

# *********************************************************************
# FUNCTION TO CONVERT SPHERICAL HARMONIC COEFFICIENTS INTO FROM A DENSITY INTERFACE INTO A GRAVITATIONAL POTENTIAL # #EVALUATED AT r, THETA, and PHI POINTS BELOW INTERFACE at radius R up to SH degree lmax
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

    vals=np.zeros((lmax,lmax))
    G=6.67e-11
    vals=np.zeros((lmax,lmax))

    for l in range(2,3):
        for m in range(l+1):
            if m==0:
                norm=clm[0,l,m]
                vals[m,l]=(delrho*(4*np.pi*G))*((norm)*(pysh.expand.spharm_lm(l,m,thet*180/np.pi,(phi*180/np.pi)-(((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/l)))*(1))*(r*(r/R)**(l-1))*(1/(2*l+1))                  
            else:
                norm=np.linalg.norm(clm[:,l,m])
                vals[m,l]=(delrho*(4*np.pi*G))*((norm)*(pysh.expand.spharm_lm(l,m,thet*180/np.pi,(phi*180/np.pi)-(((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/m)))*(1))*(r*(r/R)**(l-1))*(1/(2*l+1))
                ##include (r/R)**l-2
    return vals



