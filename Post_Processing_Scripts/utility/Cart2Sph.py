
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
	rho=np.zeros(Ncells)
	thet=np.zeros(Ncells)
	phi=np.zeros(Ncells)
	coor_arr_sph=np.zeros([Ncells,3])    

	for m in range(Ncells):	

    
		rho[m]=((Coor_arr[m,0])**2 + (Coor_arr[m,1])**2 + (Coor_arr[m,2])**2)**(1/2)
		thet[m]=np.arctan2((((Coor_arr[m,0])**2 + (Coor_arr[m,1])**2)**(1/2)),Coor_arr[m,2])
		phi[m]=np.arctan2(Coor_arr[m,1],Coor_arr[m,0])
		coor_arr_sph[m,:]=[rho[m],thet[m],phi[m]]

	return coor_arr_sph



