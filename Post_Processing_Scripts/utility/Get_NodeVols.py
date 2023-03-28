

# *********************************************************************
# FUNCTION TO EVALUATE NODE VOLUMES/AREAS FROM MESH INFORMATION
#
# *********************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd


def main(NodesData,CellsData,Cellvol):

	Nodesize_arr=len(NodesData[:,0])
	Cellsize_arr=len(CellsData[:,0])
	NodeVol=np.zeros([Nodesize_arr,1])
	Counter=np.zeros([Nodesize_arr,1])

	for y in range(Cellsize_arr): 
		for z in range(3):
			NodeVol[CellsData[y,z]]=NodeVol[CellsData[y,z]]+Cellvol[y]

	return NodeVol


