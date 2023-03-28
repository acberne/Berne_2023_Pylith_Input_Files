

# *********************************************************************
# FUNCTION TO PARSE VTK FILE WITH 2D CELLS
#
# *********************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd


def main(finame):	

	my_cols=['1' ,'2' ,'3' ,'4' ,'5', '6']

	Vtkfile=pd.read_csv(str(finame)+".vtk", 
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

	return NodesData,CellsData,DispData


