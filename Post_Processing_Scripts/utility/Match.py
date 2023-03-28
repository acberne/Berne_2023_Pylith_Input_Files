

# *********************************************************************
# FUNCTION TO SUPERIMPOSE NODES DISP or FORCE INFO FROM TWO DIFFERENT DISP OUTPUTS INTO ONE
#
# *********************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd


def main(RawData1,RawData2):	

	All_Points=RawData2[:,0:3]
	Inner_Points=RawData1[:,0:3]

	n_points=len(Inner_Points[:,0])
	SurfLoc=np.zeros(n_points)

	for i in range(n_points):
	    SurfLoc[i]=np.argmin(((All_Points[:,0]-Inner_Points[i,0])**2 + (All_Points[:,1]-Inner_Points[i,1])**2 + (All_Points[:,2]-Inner_Points[i,2])**2)**(1/2)) 
	SurfLoc=SurfLoc.astype(int)

	for i in range(n_points):
	    RawData2[SurfLoc[i],3]=RawData2[SurfLoc[i],3]+RawData1[i,3]
	    RawData2[SurfLoc[i],4]=RawData2[SurfLoc[i],4]+RawData1[i,4]
	    RawData2[SurfLoc[i],5]=RawData2[SurfLoc[i],5]+RawData1[i,5]
	    
	CombData=RawData2

	return CombData


