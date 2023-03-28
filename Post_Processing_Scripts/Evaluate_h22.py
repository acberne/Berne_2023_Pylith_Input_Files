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



ts=pd.read_csv("eccentricity.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

ecc=float(ts)


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


ts=pd.read_csv("Period.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

T=float(ts)


ts=pd.read_csv("Planet_Mass.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

Ms=float(ts)


ts=pd.read_csv("Semi_Maj.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=ts
ts=df3.to_numpy()

a=float(ts)


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


















##Step 1 Read the file containing displacement

my_cols=['hash' ,'vtk' ,'DataFile' ,'Version' ,'2.0', 'hash2']

Vtkfile=pd.read_csv("OuterSurf_AddForce.vtk", 
                          names=my_cols,
                          sep=" ",
                          index_col=False)
df1=Vtkfile
RawData=df1.to_numpy()


vtklength=len(RawData[:,1])

for x in range(vtklength):
    if RawData[x,0]=='POINTS':
        PointsStart=x+1;
        break
    
for x in range(vtklength):
    if RawData[x,0]=='CELLS':
        CellsStart=x+1;
        PointsEnd=x;
        break

for x in range(vtklength):
    if RawData[x,0]=='CELL_TYPES':
        CellsEnd=x;
        break
        
        
for x in range(vtklength):
    if RawData[x,0]=='VECTORS':
        Vectorstart=x+1;
        break


print("Evaluating h22...")


NodesData=RawData[PointsStart:PointsEnd,:]
CellsData=RawData[CellsStart:CellsEnd,:]
DispData=RawData[Vectorstart:vtklength,:]

##Reformatting
CellsData=CellsData[:,2:len(CellsData[1,:])-1]
CellsData=CellsData.astype(int)

NodesData=NodesData[:,0:3]
NodesData=NodesData.astype(float)

DispData=DispData[:,0:3]
DispData=DispData.astype(float)


def Cart2Sph(vec_arr,disp_arr):

    Ncells=len(vec_arr[:,0])

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
    RotVec=np.zeros([Ncells,3])
    vec_arr_sph=np.zeros([Ncells,3])
    
    for m in range(Ncells):
        rho[m]=((vec_arr[m,0])**2 + (vec_arr[m,1])**2 + (vec_arr[m,2])**2)**(1/2)
        thet[m]=np.arctan2((((vec_arr[m,0])**2 + (vec_arr[m,1])**2)**(1/2)),vec_arr[m,2])
        phi[m]=np.arctan2(vec_arr[m,1],vec_arr[m,0])
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
        RotVec[m,:]=np.matmul(Jac[m,:,:],disp_arr[m,:])
        ##Also want deviatoric stress
    return RotVec

DispSph=Cart2Sph(NodesData,DispData)



##Do the same with NodesData

xpoints=NodesData[:,0]
ypoints=NodesData[:,1]
zpoints=NodesData[:,2]

##Converting
veclen=len(xpoints);

Coor_arr_cart=np.zeros((veclen,3))

for i in range(veclen):
    Coor_arr_cart[i,0]=xpoints[i]
    Coor_arr_cart[i,1]=ypoints[i]
    Coor_arr_cart[i,2]=zpoints[i]


Ncells=len(Coor_arr_cart[:,0])

##Writing out conversion
rhof=np.zeros(Ncells)
thetf=np.zeros(Ncells)
phif=np.zeros(Ncells)
coor_arr_sph=np.zeros([Ncells,3])

for m in range(Ncells):
    rhof[m]=((Coor_arr_cart[m,0])**2 + (Coor_arr_cart[m,1])**2 + (Coor_arr_cart[m,2])**2)**(1/2)
    thetf[m]=np.arctan2((((Coor_arr_cart[m,0])**2 + (Coor_arr_cart[m,1])**2)**(1/2)),Coor_arr_cart[m,2])
    phif[m]=np.arctan2(Coor_arr_cart[m,1],Coor_arr_cart[m,0])
    coor_arr_sph[m,:]=[rhof[m],thetf[m],phif[m]]

    
##Want to get r thet and phi
NodesSph=coor_arr_sph



#plt.subplot(221)

#plt.figure(1)
#plt.scatter(NodesSph[:,2],NodesSph[:,1],c=DispSph[:,0])


##Now need to use scipy to interpolate data

from scipy.interpolate import griddata

#plt.subplot(222)

grid_x, grid_y = np.mgrid[0:np.pi:100j, -np.pi:np.pi:200j]

points=np.transpose([NodesSph[:,1],NodesSph[:,2]])

values=DispSph[:,0]

grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')

#plt.imshow(grid_z0.T,extent=(0,1,0,1),origin='lower')




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








grid_z0=grid_z0+OuterRad


##Pyshtools stuff


G=6.67408e-11
M=1.08e20 ##kg

MG=G*M
gm=MG
R=OuterRad ##m
 ##Test Radius Value


shape_grid=pysh.SHGrid.from_array(grid_z0)



#fig1, ax1 = shape_grid.plot(colorbar='bottom',
                           #cb_label='Topography, m')

##Define a few constants for Enceladus


shape=shape_grid.expand()

clm=shape.to_array() ##Spherical Harmonic Coefficients for gravity



def GetShape(r,thet,phi,lmax,clm):

    vals=np.zeros((lmax,lmax))

    for l in range(2,3):
        for m in range(l+1):
            if m==0:
                norm=clm[0,l,m]
                vals[m,l]=((norm)*(pysh.expand.spharm_lm(l,m,thet*180/np.pi,(phi*180/np.pi)-(((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/l)))*(1))#*(r*(r/R)**(l-1)))*(1/(2*l+1))                  
            else:
                norm=np.linalg.norm(clm[:,l,m])
                vals[m,l]=((norm)*(pysh.expand.spharm_lm(l,m,thet*180/np.pi,(phi*180/np.pi)-(((np.arctan2(clm[1,l,m],clm[0,l,m])*180/np.pi)))*(1/m)))*(1))#*(r*(r/R)**(l-1)))*(1/(2*l+1))
                ##include (r/R)**l-2


    return vals



##Evaluating Potential at grid of points



G=6.674e-11
Ms=Ms ##Mass of Saturn in kg
a=a ##semi-major axis of Enceladus about Saturn
e=ecc ## Enceladus eccentricity
T=T ##Period in seconds
ne=(np.pi * 2)/T ##Mean motion in rad/s
r=OuterRad ##Test radius, radius of Enceladus
t=(1/2)*1.370218*86400 ## Test time, 0 corresonds to periapse
prefac=-(3/2)*(G*Ms)*(1/a**3)*e*np.cos(ne*t)
veclen=30;

##Writing out acceleration

thet=np.linspace(0,np.pi,veclen)
phi=np.linspace(0,2*np.pi,veclen)
omega=2*np.pi/T
lmax=49

##This treats two coordinates, but could be extended to three if needed

Shape20=np.zeros((veclen,veclen))
Shape22=np.zeros((veclen,veclen))

for i in range(veclen):
    for j in range(veclen):

        Shape20[i,j]=GetShape(r,thet[i],phi[j],lmax,clm)[0,2]
        Shape22[i,j]=GetShape(r,thet[i],phi[j],lmax,clm)[2,2]
        





tC=pd.read_csv("Current_Time.txt", 
                          sep=" ",
                          index_col=False,
                          header=None)
df3=tC
tC=df3.to_numpy()

tC=float(tC)



G=6.674e-11
Ms=Ms ##Mass of Saturn in kg
a=a ##semi-major axis of Enceladus about Saturn
e=ecc ## Enceladus eccentricity
T=T ##Period in seconds
ne=(np.pi * 2)/T ##Mean motion in rad/s
 ##Test radius, radius of Enceladus
t=tC*T ## Test time, 0 corresonds to periapse
rho=rhoi
omega=2*np.pi/T

Vn=np.zeros((veclen,veclen))
Vr=np.zeros((veclen,veclen))

g=gt ##m/s^2

prefac=(omega**2)*(r**2)*(1/g)

##This treats two coordinates, but could be extended to three if needed

for i in range(veclen):
    for j in range(veclen):

        Vn[i,j]=(np.cos(omega*t)*(e/2)*(9/2)*np.cos(2*phi[j])*(np.sin(thet[i]))**2)  +3*e*np.sin(thet[i]) * np.sin(thet[i]) * np.sin(omega*t)* np.sin(2*phi[j])

Vr=Vn*prefac




##Laying out the forcing function



InducedV=Shape22

#h2eval=(InducedV)/(Vr)




###Grid-searching to find best value of Vr.


numtest=1000000
h2_test=np.linspace(0,2.5,numtest)


from numpy import linalg as LA


Diff=np.zeros(numtest)

for i in range (numtest):
    IndV=h2_test[i]*Vr
    Diff[i]=LA.norm(np.abs(IndV-InducedV))
    

h2=h2_test[np.argmin(Diff)]
    
    
#h2=np.median(h2eval)

file_object=open('h22_faulted.txt','a')
h2_str=str(h2)
file_object.write(h2_str)
file_object.write('\n')
file_object.close()

