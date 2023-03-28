

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

	Ncells=len(Vec_arr[:,0])

	##Writing out rotation
	Jaci=np.zeros([Ncells,3,3])
	Jac=np.zeros([Ncells,3,3])
	JacT=np.zeros([Ncells,3,3])
	rho=np.zeros(Ncells)
	thet=np.zeros(Ncells)
	phi=np.zeros(Ncells)
	xcor=np.zeros(Ncells)
	ycor=np.zeros(Ncells)
	zcor=np.zeros(Ncells)
	rhat=np.zeros(Ncells)
	thethat=np.zeros(Ncells)
	phihat=np.zeros(Ncells)
	J11=np.zeros(Ncells)
	J12=np.zeros(Ncells)
	J13=np.zeros(Ncells)
	J21=np.zeros(Ncells)
	J22=np.zeros(Ncells)
	J23=np.zeros(Ncells)
	J31=np.zeros(Ncells)
	J32=np.zeros(Ncells)
	J33=np.zeros(Ncells)
	vec_arr_cart=np.zeros([Ncells,3])
	coor_arr_cart=np.zeros([Ncells,3])
	vec_arr_sph=np.zeros([Ncells,3])
	coor_arr_sph=np.zeros([Ncells,3])

	for m in range(Ncells):
		rho[m]=Coor_arr[m,0]
		thet[m]=Coor_arr[m,1]
		phi[m]=Coor_arr[m,2]
		rhat[m]=Vec_arr[m,0]
		thethat[m]=Vec_arr[m,1]
		phihat[m]=Vec_arr[m,2]
		J11[m]=np.sin(thet[m]) * np.cos(phi[m])
		J12[m]=np.sin(thet[m]) * np.sin(phi[m])
		J13[m]=np.cos(thet[m])
		J21[m]=np.cos(thet[m])*np.cos(phi[m])
		J22[m]=np.cos(thet[m])*np.sin(phi[m])
		J23[m]=-np.sin(thet[m])
		J31[m]=-np.sin(phi[m])
		J32[m]=np.cos(phi[m])
		J33[m]=0;
		Jac[m,0,0]=J11[m]
		Jac[m,0,1]=J12[m]
		Jac[m,0,2]=J13[m]
		Jac[m,1,0]=J21[m]
		Jac[m,1,1]=J22[m]
		Jac[m,1,2]=J23[m]
		Jac[m,2,0]=J31[m]
		Jac[m,2,1]=J32[m]
		Jac[m,2,2]=J33[m] 
		Jaci[m,:,:]=np.linalg.inv(Jac[m,:,:])
		JacT[m,:,:]=np.transpose(Jac[m,:,:])
		vec_arr_sph[m,:]=[rhat[m],thethat[m],phihat[m]]
		vec_arr_cart[m,:]=np.matmul(Jaci[m,:,:],np.transpose(vec_arr_sph[m,:]))
		coor_arr_sph[m,:]=[rho[m],thet[m],phi[m]]
		xcor[m]=rho[m]*np.cos(phi[m])*np.sin(thet[m])
		ycor[m]=rho[m]*np.sin(phi[m])*np.sin(thet[m])
		zcor[m]=rho[m]*np.cos(thet[m])
		coor_arr_cart[m,:]=[xcor[m],ycor[m],zcor[m]]

	return coor_arr_cart,vec_arr_cart


