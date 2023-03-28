

# *********************************************************************
# FUNCTION TO EVALUATE CELL VOLUMES (3D) FROM MESH INFORMATION
#
# *********************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd


##Need function for evaluating determinant of tet element

def determinant_3x3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
	    m[1][0] * (m[0][1] * m[2][2] - m[0][2] * m[2][1]) +
	    m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]))


def subtract(a, b):
    return (a[0] - b[0],
	    a[1] - b[1],
	    a[2] - b[2])

def tetrahedron_calc_volume(a, b, c, d):
    return (abs(determinant_3x3((subtract(a, b),
	                         subtract(b, c),
	                         subtract(c, d),
	                         ))) / 6.0)


def main(CellsData,NodesData):


		
	Cellsize_arr=len(CellsData[:,2]);
	Cellloc=np.zeros([Cellsize_arr,3]); ##Three Sptial Dimensions
	Cellvol=np.zeros([Cellsize_arr,1]);
	Nodesize_arr=len(NodesData[:,2]);

	for x in range(Cellsize_arr):
	    Node1num= CellsData[x,0];
	    Node2num= CellsData[x,1];
	    Node3num= CellsData[x,2];
	    Node4num= CellsData[x,3];
	    Node1xval=NodesData[Node1num,0];
	    Node1yval=NodesData[Node1num,1];
	    Node1zval=NodesData[Node1num,2];
	    Node2xval=NodesData[Node2num,0];
	    Node2yval=NodesData[Node2num,1];
	    Node2zval=NodesData[Node2num,2];
	    Node3xval=NodesData[Node3num,0];
	    Node3yval=NodesData[Node3num,1];
	    Node3zval=NodesData[Node3num,2];
	    Node4xval=NodesData[Node4num,0];
	    Node4yval=NodesData[Node4num,1];
	    Node4zval=NodesData[Node4num,2];
	    xloc=(Node1xval+Node2xval+Node3xval+Node4xval)/4 ; 
	    yloc=(Node1yval+Node2yval+Node3yval+Node4yval)/4 ;
	    zloc=(Node1zval+Node2zval+Node3zval+Node4zval)/4 ;
	    ##Calculate Location
	    Cellloc[x,0]=xloc;
	    Cellloc[x,1]=yloc;
	    Cellloc[x,2]=zloc;
	    ##Calculate Volume in m3 of each tetrahedral element
	    am=[Node1xval,Node1yval,Node1zval]
	    bm=[Node2xval,Node2yval,Node2zval]
	    cm=[Node3xval,Node3yval,Node3zval]
	    dm=[Node4xval,Node4yval,Node4zval]
	    Cellvol[x]=tetrahedron_calc_volume(am, bm, cm, dm)

	    


	return Cellvol,Cellloc


