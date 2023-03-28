



# *********************************************************************
# FUNCTION TO DEFINE PARAMETERS AT START OF SIMULATION
# *********************************************************************


import numpy as np
import pandas as pd
import pickle 


###BEGIN USER INPUT SECTION

Thickness=25000 ##Shell Thickness
outrad=252000 ##Outer Radius
ecc=0.0047 ##eccentricity
rhoi=925  ##Density of Ice
rhow=1007  ##Density of Water
T=1.370218*(24*3600) ##Period in seconds
Ms=5.683e26 ##Mass of Saturn in kg
a=237.948e6 ##Semi Major Axis
gt=0.113  ##Gravity at top of ice shell
gb=0.120 ##Gravity at bottom of ice shell
Vs=1888.8  ##Ice Shear Wave Velocity
Vp=3750 ##Ice P Wave Velocity
lmax_grav=40 ##Maximum spherical harmonic degree evaluated for Eval_grav function
Size_factor=2 ##Factor to adjust mesh size
Fault_fac=5 ##Additional Factor for Scaling Faults
Static_Coefficient_Fric=0.8 #Static Coefficient of Friction (If Friction is Turned on)

##Conversion from shear and bulk to wave speeds

k=8.6e9 ## Weak Zone Bulk Modulus
G=3.3e9 ##Weak Zone Shear Modulus
WZ_Vp=((k+4*(G)/3)/(rhoi))**(1/2)
WZ_Vs=(G/rhoi)**(1/2)
Dyn_intensity=20 ##Intensity of Dynamic Remeshing
Time_inter=1##No. of timesteps



###END USER INPUT SECTION

##Define Class and Sat Object for numerical parameters and binary settings

Show_Timeleft=Show_Gui

class Sat_Template(object):
	def __init__(Satellite, outer_radius, thickness, eccentricity, rho_ice, rho_water, period, mass_primary, semi_major_axis, surface_gravity, ice_ocean_gravity,vs,vp,lmax,size_factor,fault_factor,static_coeff_fric,time_inter,wz_vs,wz_vp,dyn_intensity):

		Satellite.outer_radius = outer_radius
		Satellite.thickness = thickness
		Satellite.eccentricity = eccentricity
		Satellite.rho_ice = rho_ice
		Satellite.rho_water = rho_water
		Satellite.period = period
		Satellite.mass_primary = mass_primary
		Satellite.semi_major_axis = semi_major_axis
		Satellite.surface_gravity = surface_gravity
		Satellite.ice_ocean_gravity = ice_ocean_gravity
		Satellite.vs = vs
		Satellite.vp = vp
		Satellite.lmax = lmax
		Satellite.size_factor = size_factor
		Satellite.fault_factor = fault_factor
		Satellite.static_coeff_fric = static_coeff_fric
		Satellite.time_inter = time_inter
		Satellite.wz_vs=wz_vs
		Satellite.wz_vp=wz_vp
		Satellite.dyn_intensity=dyn_intensity


with open('Satellite.pkl', 'wb') as outp:
    Sat_Model =  Sat_Template(outrad,Thickness,ecc,rhoi,rhow,T,Ms,a,gt,gb,Vs,Vp,lmax_grav,Size_factor,Fault_fac,Static_Coefficient_Fric,Time_inter,WZ_Vs,WZ_Vp,Dyn_intensity)
    pickle.dump(Sat_Model, outp, pickle.HIGHEST_PROTOCOL)

del Sat_Model



OutRad=np.array([outrad])
inrad=outrad-Thickness
InRad=np.array([inrad])
np.savetxt("Inner_Radius.txt",InRad,fmt ="%f")
np.savetxt("Outer_Radius.txt",OutRad,fmt ="%f")
tstr=str(outrad)
print("Outer Radius:")
print(tstr)
print("m")
tstr=str(inrad)
print("Inner Radius:")
print(tstr)
print("m")





Static_Coefficient_Fric=np.array([Static_Coefficient_Fric])
np.savetxt("Static_Coefficient_Fric.txt",Static_Coefficient_Fric,fmt ="%f")
tstr=str(Static_Coefficient_Fric)
print("Static Coefficient of Friction (If Friction is Switched On):")
print(tstr)



Size_factor=np.array([Size_factor])
np.savetxt("Size_factor.txt",Size_factor,fmt ="%f")
tstr=str(Size_factor)
print("Size_factor:")
print(tstr)


Fault_fac=np.array([Fault_fac])
np.savetxt("Fault_fac.txt",Fault_fac,fmt ="%f")
tstr=str(Fault_fac)
print("Fault Scaling Factor:")
print(tstr)


Time_inter=np.array([Time_inter])
np.savetxt("Time_inter.txt",Time_inter,fmt ="%d")
tstr=str(Time_inter)
print("Time Intervals:")
print(tstr)


Dyn_intensity=np.array([Dyn_intensity])
np.savetxt("Dyn_intensity.txt",Dyn_intensity,fmt ="%d")
tstr=str(Dyn_intensity)
print("Dynamic Remeshing Intensity:")
print(tstr)

ecc=np.array([ecc])
np.savetxt("eccentricity.txt",ecc,fmt ="%f")
tstr=str(ecc)
print("Eccentricity:")
print(tstr)





Vs=np.array([Vs])
np.savetxt("Vs.txt",Vs,fmt ="%f")
tstr=str(Vs)
print("S-Wave Velocity:")
print(tstr)
print("m/s")

Vp=np.array([Vp])
np.savetxt("Vp.txt",Vp,fmt ="%f")
tstr=str(Vp)
print("P-Wave Velocity:")
print(tstr)
print("m/s")




rhoi=np.array([rhoi])
np.savetxt("Ice_Density.txt",rhoi,fmt ="%f")
tstr=str(rhoi)
print("Ice Density:")
print(rhoi)
print("kg/m^3")






rhow=np.array([rhow])
np.savetxt("Water_Density.txt",rhow,fmt ="%f")
tstr=str(rhow)
print("Water Density:")
print(rhow)
print("kg/m^3")






T=np.array([T])
np.savetxt("Period.txt",T,fmt ="%f")
T=str(T/3600)
print("Period:")
print(T)
print("Hours")





Ms=np.array([Ms])
np.savetxt("Planet_Mass.txt",Ms,fmt ="%f")
Ms=str(Ms)
print("Mass of Parent:")
print(Ms)
print("kg")



a=np.array([a])
np.savetxt("Semi_Maj.txt",a,fmt ="%f")
a=str(a)
print("Semi-Major Axis:")
print(a)
print("m")




WZ_Vs=np.array([WZ_Vs])
np.savetxt("WZ_Vs.txt",WZ_Vs,fmt ="%f")
WZ_Vs=str(WZ_Vs)
print("Weak Zone Shear Wave Speed:")
print(WZ_Vs)
print("m/s")



gt=np.array([gt])
np.savetxt("Surf_grav.txt",gt,fmt ="%f")
gt=str(gt)
print("Surface Gravity:")
print(gt)
print("m/s^2")




gb=np.array([gb])
np.savetxt("Base_grav.txt",gb,fmt ="%f")
gb=str(gb)
print("Ice-Ocean Gravity:")
print(gb)
print("m/s^2")



lmax_grav=np.array([lmax_grav])
np.savetxt("lmax_grav.txt",lmax_grav,fmt ="%f")
lmax_grav=str(lmax_grav)
print("Maximum Gravity Spherical Harmonic Degree:")
print(lmax_grav)
print("")
