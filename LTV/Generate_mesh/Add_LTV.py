##Steps: 1. Mesh a coarse-ish sphere in CUBIT. NO UNTANGLING OR SMOOTHING

## 2. Export the mesh and run it through this program

## 3. Send it back to cubit and do a finer mesh

## 4. Run Simulation!


finame = "Base_Sphere.e"

import sys
import numpy as np
import netCDF4
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
from scipy.interpolate import griddata

####################Function Definitions



def Cart2Sph(Coor_arr):

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


def Sph2Cart(Coor_arr):	
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




def Interp_Grid(NodesSph):

	grid_x, grid_y = np.mgrid[0:np.pi:300j, -np.pi:np.pi:300j]

	points=np.transpose([NodesSph[:,1],NodesSph[:,2]])

	values=NodesSph[:,0]

	grid_z0 = griddata(points, values, (grid_x, grid_y), method='cubic',fill_value=-1)

	grid_fill = griddata(points, values, (grid_x, grid_y), method='nearest')


	for i in range(grid_z0.shape[0]):
		for j in range(grid_z0.shape[1]):
			if grid_z0[i,j]==-1:
				grid_z0[i,j]=grid_fill[i,j]



	return grid_z0


def clm2Topo_comp(theta,phi,lmax,clm,fac):
	n_co=1.0
	size=np.size(phi)
	vals=np.zeros((lmax,lmax,size))
	vec_arr=np.array([phi,theta])
	def GetTopo(theta,phi,lmax,clm):
		for l in range(2,lmax):
			for m in range(l+1):
				if m==0:
					norm=clm[0,l,m]*(1/(l-1)**n_co)
					vals[m,l,:]=((norm*(pysh.expand.spharm_lm(l,m,theta*180/np.pi,phi*180/np.pi-((((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/l))))))
				else:
					norm=np.linalg.norm(clm[:,l,m])*(1/(l-1)**n_co)
					vals[m,l,:]=((norm*(pysh.expand.spharm_lm(l,m,theta*180/np.pi,phi*180/np.pi-((((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/m))))))
				##include (r/R)**l-2
		return vals
	WorksDown=GetTopo(theta,phi,lmax,clm)
	sums=np.sum(WorksDown,axis=1)
	sums2=np.sum(sums,axis=0)
	Topovals=sums2
	Topovals=Topovals*fac
	return Topovals

def clm2Topo(theta,phi,lmax,clm,fac):

	size=np.size(phi)
	vals=np.zeros((lmax,lmax,size))
	vec_arr=np.array([phi,theta])
	def GetTopo(theta,phi,lmax,clm):
		for l in range(2,lmax):
			for m in range(l+1):
				if m==0:
					norm=clm[0,l,m]
					vals[m,l,:]=((norm*(pysh.expand.spharm_lm(l,m,theta*180/np.pi,phi*180/np.pi-((((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/l))))))
				else:
					norm=np.linalg.norm(clm[:,l,m])
					vals[m,l,:]=((norm*(pysh.expand.spharm_lm(l,m,theta*180/np.pi,phi*180/np.pi-((((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/m))))))
				##include (r/R)**l-2
		return vals
	WorksDown=GetTopo(theta,phi,lmax,clm)
	sums=np.sum(WorksDown,axis=1)
	sums2=np.sum(sums,axis=0)
	Topovals=sums2
	Topovals=Topovals*fac
	return Topovals

#######################################
##Import Params


ts=pd.read_csv("fact.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

fact=float(ts)

ts=pd.read_csv("Ice_Density.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

rhoi=float(ts)

ts=pd.read_csv("Water_Density.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

rhow=float(ts)



ts=pd.read_csv("Surf_grav.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

gt=float(ts)


ts=pd.read_csv("Base_grav.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

gb=float(ts)





######################################
##First, import the mesh

# ----------------------------------------------------------------------
# Get coordinates of points from ExodusII file.
exodus = netCDF4.Dataset(finame, 'a')
coords = np.empty((len(exodus.variables['coordx'][:]),3))
coords[:,0] = exodus.variables['coordx'][:]
coords[:,1] = exodus.variables['coordy'][:]
coords[:,2] = exodus.variables['coordz'][:]

# ----------------------------------------------------------------------
# Rotate into spherical


NodesData=coords
NodesSph=Cart2Sph(NodesData)



###########################################################



def readtxt(finame):
        my_cols=['hash' ,'vtk' ,'DataFile' ,'Version','Another','Column','last']
        filevar=pd.read_csv(str(finame)+".txt", 
                                sep=",",
                                index_col=False,
                                header=None)
        df1=filevar
        var=df1.to_numpy()
        return var



LFT=readtxt('Edited')
np.savetxt("uh",LFT)
LFT=LFT*1e3 ##Convert to m
LFT=LFT.astype(float)
ShapeSph=Cart2Sph(LFT) 
Shape_grid=Interp_Grid(ShapeSph)
#plt.imshow(Shape_grid)
#plt.colorbar()
#plt.show()



###################################################################




ts=pd.read_csv("Outer_Radius.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()
OutRad=float(ts)


#shape_grid=pysh.SHGrid.from_ellipsoid(lmax=100,a=arad+fac,b=brad+fac,c=crad-fac,grid='GLQ')
#shape=shape_grid.expand()
#clm_ref=shape.to_array() 
#Ref_Topo=clm2Topo(ShapeSph[:,1],ShapeSph[:,2],3,clm_ref)
#Ref_Ellipsoid=ShapeSph
#Ref_Ellipsoid[:,0]=Ref_Topo+OutRad

#ShapeSph=Cart2Sph(LFT) 
#Residual[:,0]=ShapeSph[:,0]-Residual[:,0]
#Residual_grid=Interp_Grid(Residual)

#plt.imshow(Residual_grid)
#plt.colorbar()
#plt.show()



###########################Let's go a different Route


ShapeSph=Cart2Sph(LFT)
shape_grid=pysh.SHGrid.from_array(Shape_grid)
shape=shape_grid.expand()
clm=shape.to_array() 
clm_nhyd=clm
clm_nhyd[0,2,0]= clm_nhyd[0,2,0]*(1-((2.569)/3.846))
clm_nhyd[0,2,2]= clm_nhyd[0,2,2]*(1-((771)/917)) ##Subtracting out the hydrostatic components 
#clm_nhyd[0,0,0]= 0
fac=1
nhyd_Topo=clm2Topo_comp(NodesSph[:,1],NodesSph[:,2],30,clm_nhyd,fac)
nhyd_Shape=NodesSph
nhyd_Shape[:,0]=nhyd_Topo
nhyd_Shape_grid=Interp_Grid(nhyd_Shape)
#plt.imshow(nhyd_Shape_grid)
#plt.colorbar()
#plt.show()


##Long Wavelengths

ShapeSph=Cart2Sph(LFT)
shape_grid=pysh.SHGrid.from_array(Shape_grid)
shape=shape_grid.expand()
clm=shape.to_array() 
clm_nhyd=clm
clm_nhyd[0,2,0]= clm_nhyd[0,2,0]*(1-((2.569)/3.846))
clm_nhyd[0,2,2]= clm_nhyd[0,2,2]*(1-((771)/917)) ##Subtracting out the hydrostatic components 
#clm_nhyd[0,0,0]= 0
nhyd_Topo_2=clm2Topo(NodesSph[:,1],NodesSph[:,2],5,clm_nhyd,0)
nhyd_Shape=NodesSph
nhyd_Shape[:,0]=nhyd_Topo
nhyd_Shape_grid=Interp_Grid(nhyd_Shape)
#plt.imshow(nhyd_Shape_grid)
#plt.colorbar()
#plt.show()



##Top Surface


ShapeSph=Cart2Sph(LFT)
shape_grid=pysh.SHGrid.from_array(Shape_grid)
shape=shape_grid.expand()
clm=shape.to_array() 
clm_nhyd=clm
clm_nhyd[0,2,0]= clm_nhyd[0,2,0]*(1-((2.569)/3.846))
clm_nhyd[0,2,2]= clm_nhyd[0,2,2]*(1-((771)/917)) ##Subtracting out the hydrostatic components 
#clm_nhyd[0,0,0]= 0
nhyd_Topo_3=clm2Topo(NodesSph[:,1],NodesSph[:,2],30,clm_nhyd,1)
nhyd_Shape=NodesSph
nhyd_Shape[:,0]=nhyd_Topo
nhyd_Shape_grid=Interp_Grid(nhyd_Shape)
#plt.imshow(nhyd_Shape_grid)
#plt.colorbar()
#plt.show()





# ----------------------------------------------------------------------
# Modify Inner Surface

rhow=rhow

rhoi=rhoi


ts=pd.read_csv("Outer_Radius.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

OuterRad=float(ts)



ts=pd.read_csv("Inner_Radius.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

InnerRad=float(ts)
Thic=OuterRad-InnerRad
Topovals=nhyd_Topo+nhyd_Topo_2

##Getting Hydrostatic component


ShapeSph=Cart2Sph(LFT)
shape_grid=pysh.SHGrid.from_array(Shape_grid)
shape=shape_grid.expand()
clm=shape.to_array() 
shape_grid=pysh.SHGrid.from_ellipsoid(lmax=100,a=OuterRad,b=OuterRad,c=OuterRad,grid='GLQ')
shape=shape_grid.expand()
clm_hyd=shape.to_array() 
clm_hyd[0,2,0]= clm[0,2,0]*(((2.569)/3.846))
clm_hyd[0,2,2]= clm[0,2,2]*(((771)/917)) ##Subtracting out the hydrostatic components 
#clm_hyd[0,0,0]= 0
hyd_Topo=clm2Topo(NodesSph[:,1],NodesSph[:,2],5,clm_hyd,1)
hyd_Shape=NodesSph
hyd_Shape[:,0]=hyd_Topo
hyd_vals=hyd_Topo


NodesSph=Cart2Sph(NodesData)

#print(NodesSph[:,0])

#print(Topovals)
Checkmat=np.zeros(len(NodesSph[:,0]))
LocTop=np.zeros(len(NodesSph[:,0]))
LocBottom=np.zeros(len(NodesSph[:,0]))
OldLocTop=np.zeros(len(NodesSph[:,0]))
OldLocBottom=np.zeros(len(NodesSph[:,0]))

for i in range(len(NodesSph[:,0])):
    if NodesSph[i,0]<(InnerRad+300):
        NodesSph[i,0]=NodesSph[i,0]-((Topovals[i])*(rhoi/(rhow-rhoi)))*(gt/gb)*((OuterRad**2)/(InnerRad**2))*1.0+hyd_vals[i]  ##Airy Isostasy
    elif NodesSph[i,0]>(OuterRad-500):
        NodesSph[i,0]=NodesSph[i,0]+nhyd_Topo_3[i]+hyd_vals[i]  
    else:
        NodesSph[i,0]=NodesSph[i,0]
        Checkmat[i]=-1 ##Flag

for i in range(len(Checkmat)):
	if Checkmat[i] == -1:
		LocTop[i]=nhyd_Topo_3[i]+hyd_vals[i]+OuterRad#-NodesSph[i,0] 
		LocBottom[i]=  InnerRad-((Topovals[i])*(rhoi/(rhow-rhoi)))*(gt/gb)*((OuterRad**2)/(InnerRad**2))*1.0+hyd_vals[i]
		#OldLocTop[i]=(OuterRad-NodesSph[i,0])/((OuterRad-InnerRad)) #nhyd_Topo_3[i]+hyd_vals[i]+OuterRad#-NodesSph[i,0] 
		OldLocBottom[i]= (NodesSph[i,0]-InnerRad)/((OuterRad-InnerRad))#-((Topovals[i])*(rhoi/(rhow-rhoi)))*(gt/gb)*((OuterRad**2)/(InnerRad**2))*1.0+hyd_vals[i]

		NodesSph[i,0]=(LocTop[i]-LocBottom[i])*(OldLocBottom[i]) + LocBottom[i]

		
#plt.hist(NodesSph[:,0])
#plt.show()






##Need to add back on the hydrostatic component

#for i in range(len(NodesSph[:,0])):
#	NodesSph[i,0]=NodesSph[i,0]+Topovals[i]  






# Convert Back to Cartesian


Ncells=len(NodesSph[:,0])

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
    rho[m]=NodesSph[m,0]
    thet[m]=NodesSph[m,1]
    phi[m]=NodesSph[m,2]
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
    xcor[m]=rho[m]*np.cos(phi[m])*np.sin(thet[m])
    ycor[m]=rho[m]*np.sin(phi[m])*np.sin(thet[m])
    zcor[m]=rho[m]*np.cos(thet[m])
    coor_arr_cart[m,:]=[xcor[m],ycor[m],zcor[m]]

xcoordsnew=coor_arr_cart[:,0]
ycoordsnew=coor_arr_cart[:,1]
zcoordsnew=coor_arr_cart[:,2]

##Put it all back into the exodus file

exodus.variables['coordx'][:]=xcoordsnew
exodus.variables['coordy'][:]=ycoordsnew
exodus.variables['coordz'][:]=zcoordsnew

exodus.close()
