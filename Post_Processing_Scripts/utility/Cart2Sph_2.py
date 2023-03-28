
# *********************************************************************
# FUNCTION TO CONVERT COORDINATES AND VECTORS FROM SPHERICAL TO CARTESIAN FRAME
#
# *********************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd



def main(Coor_arr,Vec_arr):

	Ncells=len(Coor_arr[:,0])

	##Writing out rotation
	Jaci=np.zeros([Ncells,3,3])
	Jac=np.zeros([Ncells,3,3])
	JacT=np.zeros([Ncells,3,3])
	Cart_Stressmat=np.zeros([Ncells,3,3])
	Sph_mat_right=np.zeros([Ncells,3,3])
	Sph_Stressmat=np.zeros([Ncells,3,3])
	Exam_Sph=np.zeros([Ncells,6])
	Dev_Stress=np.zeros([Ncells,6])
	rho=np.zeros(Ncells)
	thet=np.zeros(Ncells)
	phi=np.zeros(Ncells)
	J11=np.zeros(Ncells)
	J12=np.zeros(Ncells)
	J13=np.zeros(Ncells)
	J21=np.zeros(Ncells)
	J22=np.zeros(Ncells)
	J23=np.zeros(Ncells)
	J31=np.zeros(Ncells)
	J32=np.zeros(Ncells)
	J33=np.zeros(Ncells)
	vec_arr_sph=np.zeros([Ncells,3])
	coor_arr_sph=np.zeros([Ncells,3])    
	for m in range(Ncells):	
		rho[m]=((Coor_arr[m,0])**2 + (Coor_arr[m,1])**2 + (Coor_arr[m,2])**2)**(1/2)
		thet[m]=np.arctan2((((Coor_arr[m,0])**2 + (Coor_arr[m,1])**2)**(1/2)),Coor_arr[m,2])
		phi[m]=np.arctan2(Coor_arr[m,1],Coor_arr[m,0])
		J11[m]=np.sin(thet[m]) * np.cos(phi[m])
		J12[m]=np.sin(thet[m]) * np.sin(phi[m])
		J13[m]=np.cos(thet[m])
		J21[m]=np.cos(thet[m])*np.cos(phi[m])
		J22[m]=np.cos(thet[m])*np.sin(phi[m])
		J23[m]=-np.sin(thet[m])
		J31[m]=-np.sin(phi[m])
		J32[m]=np.cos(phi[m])
		J33[m]=0;
		Jaci[m,0,0]=J11[m]
		Jaci[m,0,1]=J12[m]
		Jaci[m,0,2]=J13[m]
		Jaci[m,1,0]=J21[m]
		Jaci[m,1,1]=J22[m]
		Jaci[m,1,2]=J23[m]
		Jaci[m,2,0]=J31[m]
		Jaci[m,2,1]=J32[m]
		Jaci[m,2,2]=J33[m] 
		Jac[m,:,:]=Jaci[m,:,:]
		vec_arr_sph[m,:]=np.matmul(Jac[m,:,:],Vec_arr[m,:])
		coor_arr_sph[m,:]=[rho[m],thet[m],phi[m]]

	return coor_arr_sph,vec_arr_sph




