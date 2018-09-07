import yt
import numpy as np
from numpy import pi
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
import time 
yt.enable_parallelism()


# ===============================================
#	Decomposition of Psi4 into spinweighted spherical harmonics  
#           		l = 2,3,4 
# 	Input: Any dataset with Weyl1, Weyl2 and Radius of sphere r
#       Output: \psi_4 *r 
# ===============================================

# ===============================================
# Radius of sphere (center is set automatically in as domain_right_edge/2)
Rad = 60
# ===============================================

# Data file
loc = '../../GRChombo/GRChombo_1.0/S1x00*'


#Loading dataset

ds = yt.load(loc)

# =================
# Init Outputarrays
# =================
timedata = []
Cycletimedata = []
Spherdata = []

time0 = ds[0].current_time
time1 = ds[1].current_time
DeltaT = float(time1-time0)


# Definitions for Quadrature sceme 

N = 89

coord = np.loadtxt('PointDistFiles/lebedev/lebedev_%03d.txt'% N)
theta = coord[:,1]*pi/180; phi = coord[:,0]*pi/180 + pi
w = coord[:,2]

Starttime = 380
Startindex = (Starttime+Rad)/DeltaT
counter = 0 

Int = np.zeros(len(coord),dtype = "complex128")

# ==============================
#         Loop over all frames 
# ==============================
for i in ds:

	if counter > Startindex : 
		startCycle = time.time()

		i.print_stats() 
		center =  (i.domain_right_edge/2)

		start = time.time()

# ==================================================
	# Initalising  


		start = time.time()
		
		index = 0 

		for (k,x) in enumerate(phi):
			phi_var = phi[k]
			theta_var = theta[k]
			x1 = Rad*np.cos(phi_var)*np.sin(theta_var)+float(center[0])
			y1 = Rad*np.sin(phi_var)*np.sin(theta_var)+float(center[1])
			z1 = Rad*np.cos(theta_var)+float(center[2])
			c = [x1,y1,z1]
			ptn = i.point(c)
			ReWeyl = float(ptn["Weyl1"][0])
			ImWeyl = float(ptn["Weyl2"][0])
			Weyl4 = ReWeyl + 1j*ImWeyl
		
			Int[k] += Weyl4*DeltaT
# ==================================================
# DATA WRITEOUT
# ==================================================

		Spher = 0 
		for (k,x) in enumerate(phi):
			Spher += 4*pi*w[k]*np.absolute(Int[k])**2
		np.savetxt('time.out',timedata)
		np.savetxt('time.out',timedata)

		Spher = Spher*(Rad)**2/(16*np.pi)
	
		# time
		timedata.append(i.current_time)
		Spherdata.append(Spher)

		np.savetxt('time.out',timedata)
		np.savetxt('Spher.out',Spherdata)

# ==================================================
# Performance 
# ==================================================
		Cycletime = abs(time.time()-startCycle)
		Cycletimedata.append(Cycletime)
		np.savetxt('Cycletime.out',Cycletimedata)
# ==================================================

	print (" ")
	print (" ")
	print (" ")
	print ("current time = ", (float(i.current_time - time0) -60))
	print ("Starttime = ",Starttime )
	print (" ")
	print (" ")
	print (" ")
	counter += 1 
print("FINISHED !!! ")


