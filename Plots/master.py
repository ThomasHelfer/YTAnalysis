import yt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
import time
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

yt.enable_parallelism()

# Plot parameters
line = 2.5 # Line width
alp = 0.6 # alpha

#Loading dataset
ds = yt.load('../*.hdf5')

def checkneg(a):
        for i in np.arange(len(a)):
                if a[i]<0 :
                        return i

        return 0

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])

def _V(field, data):
	return 0.5*data["phi"]**2

def _Kinetic(field, data):
	return 0.5*data["Pi"]**2

def _Gradient(field, data):
	return 0.5*data["phi_gradient_magnitude"]**2

def _H2(field, data):
	return data["Ham"]*data["Ham"]

Gradient_data = []
Potential_data = []
Kinetic_data = []
Rho_data = []
Efolds_data = []
Gradient_data = []
L2H_data = []
time_data = []
CycleData = []

counter = 0

for i in ds:

	counter += 1

        start = time.time()

	i.add_gradient_fields(('chombo', u'phi'))
	i.add_field("V",_V, units = "")
	i.add_field("Kinetic",_Kinetic, units = "")
	i.add_field("Gradient",_Gradient, units = "cm**(-2)")
	i.add_field("H2",_H2, units = "")

	for j in sorted(i.field_list):
  		print(j)

	ad = i.all_data()

	Rho_max = ad["rho"].max()

	time_data.append(i.current_time)


#==============================================
# 	Mean potential Energydensity
#==============================================


	Debug = time.time()

        Pot = ad.mean("V", weight="cell_volume")
	Potential_data.append(Pot)

	print ("===============================================")
	print ("Time to compute Potential Energy")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("computed average value = ", Pot )
	print ("===============================================")

#==============================================
# 	Mean potential energydensity
#==============================================

	Debug = time.time()

        Kin = float(ad.mean("Kinetic", weight="cell_volume"))
	Kinetic_data.append(Kin)

	print ("===============================================")
	print ("Time to compute Kinetic Energy")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("computed average value = ", Kin )
	print ("===============================================")


#==============================================
# 	Mean gradient energydensity
#==============================================

	Debug = time.time()

        Grad = float(ad.mean("Gradient", weight="cell_volume"))
	Gradient_data.append(Grad)

	print ("===============================================")
	print ("Time to compute Kinetic Energy")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("computed average value = ", Grad )
	print ("===============================================")

#==============================================
# 	Mean Total Energy
#==============================================

	Debug = time.time()

        Rho = float(ad.mean("rho", weight="cell_volume"))
	Rho_data.append(Rho)

	print ("===============================================")
	print ("Time to compute Kinetic Energy")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("computed average value = ", Rho )
	print ("===============================================")

#==============================================
# 	Quality check
#==============================================


	print ("===============================================")
	print ("Quality check. Q = (<Kin>+<Grad>+<Pot>-<rho>)/<rho>")
	print ("Q = ", (Grad+Kin+Pot-Rho)/Rho*100, "%" )
	print ("===============================================")

#==============================================
# 	Mean Total L2H
#==============================================


	Debug = time.time()

        meanH2 = ad.mean("H2", weight="cell_volume")
	L2H = np.sqrt(meanH2)
	L2H_data.append(L2H)


	print ("===============================================")
	print ("Time to compute L2H")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("computed average value = ", L2H )
	print ("===============================================")

#==============================================
# 	Mean Total Efolds
#==============================================

	Debug = time.time()

        Efolds = float(ad.mean("Efolds", weight="cell_volume"))
	Efolds_data.append(Efolds)

	print ("===============================================")
	print ("Time to compute Efolds")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("computed average value = ", Efolds )
	print ("===============================================")

#==============================================
#	write out
#==============================================


 	CycleData.append(time.time()-start)

	np.savetxt('time.out',time_data)
	np.savetxt('Potential.out',Potential_data)
	np.savetxt('Kinetic.out',Kinetic_data)
	np.savetxt('Gradient.out',Gradient_data)
	np.savetxt('Efolds.out',Efolds_data)
	np.savetxt('rho.out',Rho_data)
	np.savetxt('L2H.out',L2H_data)
