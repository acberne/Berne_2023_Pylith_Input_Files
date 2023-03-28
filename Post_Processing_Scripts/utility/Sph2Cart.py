

# *********************************************************************
# FUNCTION TO CONVERT COORDINATES FROM SPHERICAL TO CARTESIAN FRAME
#
# *********************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd


def main(Coor_arr):	
	Ncells=len(Coor_arr[:,0])
	coor_arr_cart=np.zeros([Ncells,3])
	xcor=np.zeros(Ncells)
	ycor=np.zeros(Ncells)
	zcor=np.zeros(Ncells)
	rho=np.zeros(Ncells)
	thet=np.zeros(Ncells)
	phi=np.zeros(Ncells)
	for m in range(Ncells):	
		rho[m]=Coor_arr[m,0]
		thet[m]=Coor_arr[m,1]
		phi[m]=Coor_arr[m,2]
		xcor[m]=rho[m]*np.cos(phi[m])*np.sin(thet[m])
		ycor[m]=rho[m]*np.sin(phi[m])*np.sin(thet[m])
		zcor[m]=rho[m]*np.cos(thet[m])
		coor_arr_cart[m,:]=[xcor[m],ycor[m],zcor[m]]

	return coor_arr_cart


