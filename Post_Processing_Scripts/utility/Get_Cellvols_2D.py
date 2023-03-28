

# *********************************************************************
# FUNCTION TO EVALUATE CELL AREAS (2D) FROM MESH INFORMATION
#
# *********************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd


def det(a):
	return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[0][0]*a[1][2]*a[2][1]

	#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = det([[1,a[1],a[2]],
	     [1,b[1],b[2]],
	     [1,c[1],c[2]]])
    y = det([[a[0],1,a[2]],
	     [b[0],1,b[2]],
	     [c[0],1,c[2]]])
    z = det([[a[0],a[1],1],
	     [b[0],b[1],1],
	     [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#dot product of vectors a and b
def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

#cross product of vectors a and b
def cross(a, b):
    x = a[1] * b[2] - a[2] * b[1]
    y = a[2] * b[0] - a[0] * b[2]
    z = a[0] * b[1] - a[1] * b[0]
    return (x, y, z)	

#area of polygon poly
def area(poly):
	if len(poly) < 3: # not a plane - no area
		return 0

	total = [0, 0, 0]
	for i in range(len(poly)):
		vi1 = poly[i]
		if i is len(poly)-1:
			vi2 = poly[0]
		else:
			vi2 = poly[i+1]
		prod = cross(vi1, vi2)
		total[0] += prod[0]
		total[1] += prod[1]
		total[2] += prod[2]
	result = dot(total, unit_normal(poly[0], poly[1], poly[2]))
	return abs(result/2)




def main(CellsData,NodesData):


	

	Cellsize_arr=len(CellsData[:,2]);
	Cellloc=np.zeros([Cellsize_arr,3]); ##Three Sptial Dimensions
	Cellvol=np.zeros([Cellsize_arr,1]);
	Nodesize_arr=len(NodesData[:,2]);
	Cellvol2=np.zeros([Cellsize_arr,1]);

	for x in range(Cellsize_arr):
		Node1num= CellsData[x,0];
		Node2num= CellsData[x,1];
		Node3num= CellsData[x,2];
		Node1xval=NodesData[Node1num,0];
		Node1yval=NodesData[Node1num,1];
		Node1zval=NodesData[Node1num,2];
		Node2xval=NodesData[Node2num,0];
		Node2yval=NodesData[Node2num,1];
		Node2zval=NodesData[Node2num,2];
		Node3xval=NodesData[Node3num,0];
		Node3yval=NodesData[Node3num,1];
		Node3zval=NodesData[Node3num,2];
		xloc=(Node1xval+Node2xval+Node3xval)/3 ; 
		yloc=(Node1yval+Node2yval+Node3yval)/3 ;
		zloc=(Node1zval+Node2zval+Node3zval)/3 ;
		##Calculate Location
		Cellloc[x,0]=xloc;
		Cellloc[x,1]=yloc;
		Cellloc[x,2]=zloc;
		##Calculate Area in m^2 of each triangle
		poly = [[Node1xval, Node1yval, Node1zval], [Node2xval, Node2yval, Node2zval], [Node3xval, Node3yval, Node3zval]]
		Cellvol[x]=area(poly);


	return Cellvol,Cellloc


