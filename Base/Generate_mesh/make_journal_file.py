

# *********************************************************************
# FUNCTION TO GENERATE JOU FILE W/FAULTS FROM TXT FILE
#
# *********************************************************************

import numpy as np
import pandas as pd
from utility import Sph2Cart_2
from utility import Cart2Sph_2
from utility import Sph2Cart
from utility import Cart2Sph
from utility import Parse_txt
import pickle


class Bin_Template(object):
	def __init__():
		print() ##Null

with open('Bin_object.pkl', 'rb') as inp:
	Bin_Model = pickle.load(inp)


class Sat_Template(object):
	def __init__():
		print() ##Null

with open('Satellite.pkl', 'rb') as inp:
	Sat_Model = pickle.load(inp)
	OuterRad=Sat_Model.outer_radius
	Thickness=Sat_Model.thickness
	size_factor=Sat_Model.size_factor
	fault_factor=Sat_Model.fault_factor

InnerRad=OuterRad-Thickness

Spatialdb=open("Template.jou","w")

Spatialdb.write("reset	\n")
Spatialdb.write("set journal on	\n")
Spatialdb.write("${Outer=")
Spatialdb.write(str(OuterRad))
Spatialdb.write("}\n")
Spatialdb.write("${Inner=")
Spatialdb.write(str(InnerRad))
Spatialdb.write("}\n")
Spatialdb.write("${size=")
Spatialdb.write(str(0.001*Thickness/size_factor))
Spatialdb.write("}\n")
Spatialdb.write("${FaultWidth=")
Spatialdb.write(str(Thickness/fault_factor))
Spatialdb.write("}\n")
Spatialdb.write("create sphere radius {")
Spatialdb.write("Outer")
Spatialdb.write("} inner radius {")
Spatialdb.write("Inner")
Spatialdb.write("}\n")

##Now, we initiate a loop to write out curves, first need to import coordinates


No_faults=4

for j in range(No_faults):

	Fault_no=int(j+1)
	Data=Parse_txt.main("Fault_"+str(Fault_no))
	lat=Data[:,0]
	lon=Data[:,1]

	for i in range(len(lat)):
		if i ==0:
			Spatialdb.write("create vertex x ")
			Spatialdb.write(str((OuterRad+20000)*np.sin(((90-lat[i])*(np.pi/180)))*np.cos((360+lon[i])*(np.pi/180))))
			Spatialdb.write("	y	")
			Spatialdb.write(str((OuterRad+20000)*np.sin(((90-lat[i])*(np.pi/180)))*np.sin((360+lon[i])*(np.pi/180))))
			Spatialdb.write("	z	")
			Spatialdb.write(str((OuterRad+20000)*np.cos(((90-lat[i])*(np.pi/180)))))
			Spatialdb.write("\n")
			Spatialdb.write("${pBegin=Id('vertex')}\n")
		else:
			Spatialdb.write("create vertex x ")
			Spatialdb.write(str((OuterRad+20000)*np.sin(((90-lat[i])*(np.pi/180)))*np.cos((360+lon[i])*(np.pi/180))))
			Spatialdb.write("	y	")
			Spatialdb.write(str((OuterRad+20000)*np.sin(((90-lat[i])*(np.pi/180)))*np.sin((360+lon[i])*(np.pi/180))))
			Spatialdb.write("	z	")
			Spatialdb.write(str((OuterRad+20000)*np.cos(((90-lat[i])*(np.pi/180)))))
			Spatialdb.write("\n")

	Spatialdb.write("${pEnd=Id('vertex')}\n")
	Spatialdb.write("create curve spline vertex {pBegin} to {pEnd} delete\n")
	Spatialdb.write("${curve1=Id('curve')}\n")


	for i in range(len(lat)):
		if i ==0:
			Spatialdb.write("create vertex x ")
			Spatialdb.write(str((InnerRad-20000)*np.sin(((90-lat[i])*(np.pi/180)))*np.cos((360+lon[i])*(np.pi/180))))
			Spatialdb.write("	y	")
			Spatialdb.write(str((InnerRad-20000)*np.sin(((90-lat[i])*(np.pi/180)))*np.sin((360+lon[i])*(np.pi/180))))
			Spatialdb.write("	z	")
			Spatialdb.write(str((InnerRad-20000)*np.cos(((90-lat[i])*(np.pi/180)))))
			Spatialdb.write("\n")
			Spatialdb.write("${pBegin=Id('vertex')}\n")
		else:
			Spatialdb.write("create vertex x ")
			Spatialdb.write(str((InnerRad-20000)*np.sin(((90-lat[i])*(np.pi/180)))*np.cos((360+lon[i])*(np.pi/180))))
			Spatialdb.write("	y	")
			Spatialdb.write(str((InnerRad-20000)*np.sin(((90-lat[i])*(np.pi/180)))*np.sin((360+lon[i])*(np.pi/180))))
			Spatialdb.write("	z	")
			Spatialdb.write(str((InnerRad-20000)*np.cos(((90-lat[i])*(np.pi/180)))))
			Spatialdb.write("\n")


	Spatialdb.write("${pEnd=Id('vertex')}\n")
	Spatialdb.write("create curve spline vertex {pBegin} to {pEnd} delete\n")
	Spatialdb.write("${curve2=Id('curve')}	\n")
	Spatialdb.write("create surface skin curve {curve1} {curve2}	\n")
	Spatialdb.write("${surface=Id('surface')}	\n")
	Spatialdb.write("sweep surface {surface} vector 0.5 0.5 0.0 distance ")
	Spatialdb.write("{FaultWidth}\n")
	Spatialdb.write("${volume2=Id('volume')}	\n")
	Spatialdb.write("chop volume ")
	Spatialdb.write(str(3*j+1))
	Spatialdb.write(" with volume ")
	Spatialdb.write(str(3*j+2))
	Spatialdb.write("\n")
	Spatialdb.write("delete curve {curve1} {curve2}	\n")
	Spatialdb.write("merge volume all	\n")


Spatialdb.write("{Units('si')} 	\n")



#Uniform resolution tetmesh.



Spatialdb.write("set tetmesher HPC off	\n")
Spatialdb.write("Trimesher geometry sizing off	\n")
Spatialdb.write("surface all scheme trimesh geometry approximation angle 4 minimum size {size*km}	\n")
Spatialdb.write("volume all scheme tetmesh	\n")

Spatialdb.write(" volume 13 size {size*km}\n")
Spatialdb.write(" volume 3 6 9 12 size {FaultWidth*0.5}\n")
Spatialdb.write(" mesh volume 3 6 9 12\n")
Spatialdb.write(" mesh surface all\n")
Spatialdb.write(" mesh volume 13\n")

Spatialdb.write("volume all size {size*km}	\n")
Spatialdb.write("mesh surface all	\n")
Spatialdb.write("mesh volume all	\n")

##Need to fix this part eventually


Spatialdb.write("nodeset 1 add surface 12	\n")
Spatialdb.write("nodeset 1 name \"Fault_1\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 2 add surface 28	\n")
Spatialdb.write("nodeset 2 name \"Fault_2\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 3 add node all	\n")
Spatialdb.write("nodeset 3 name \"All_nodes\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 4 add surface 14 32 73 49 68	\n")
Spatialdb.write("nodeset 4 name \"Outer_Surf\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 5 add surface 13 31 50 67 74	\n")
Spatialdb.write("nodeset 5 name \"Inner_Surf\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 6 add curve 18 23	\n")
Spatialdb.write("nodeset 6 name \"Fault_1_edge\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 7 add curve 54 58	\n")
Spatialdb.write("nodeset 7 name \"Fault_2_edge\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 8 add curve 94 99	\n")
Spatialdb.write("nodeset 8 name \"Fault_3_edge\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 9 add curve 130 134	\n")
Spatialdb.write("nodeset 9 name \"Fault_4_edge\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 10 add surface 48	\n")
Spatialdb.write("nodeset 10 name \"Fault_3\"	\n")
Spatialdb.write("		\n")
Spatialdb.write("nodeset 11 add surface 64	\n")
Spatialdb.write("nodeset 11 name \"Fault_4\"	\n")
Spatialdb.write("export mesh \"/home/aberne/CUBIT/Coreform-Cubit-2021.11/bin/Basic_Sphere_base.exo\"  overwrite	\n")
Spatialdb.write("export mesh \"/home/aberne/CUBIT/Coreform-Cubit-2021.11/bin/Basic_Sphere_Dyn.exo\"  overwrite	\n")
Spatialdb.write("		\n")
Spatialdb.write("block 13 volume all	\n")
Spatialdb.write("export mesh \"/home/aberne/CUBIT/Coreform-Cubit-2021.11/bin/sizing.exo\" dimension 3 block 13 overwrite	\n")
















