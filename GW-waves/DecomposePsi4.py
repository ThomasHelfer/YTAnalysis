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
Rad = 120
# ===============================================

# Data file
loc = '../../GRChombo/GRChombo_1.0/S1x00*'


#Loading dataset

ds = yt.load(loc)

# =================
# Init Outputarrays
# =================
# l = 2  
Weyl4_l2_m0_data  = []
# positive m  
Weyl4_l2_m1_data = []
Weyl4_l2_m2_data = []
# negative m  
Weyl4_l2_m1n_data = []
Weyl4_l2_m2n_data = []
# l = 3  
Weyl4_l3_m0_data = []
# positive m  
Weyl4_l3_m1_data = []  
Weyl4_l3_m2_data = []
Weyl4_l3_m3_data = []
# negative m  
Weyl4_l3_m1n_data = [] 
Weyl4_l3_m2n_data = []
Weyl4_l3_m3n_data = []
# l = 4  
Weyl4_l4_m0_data = []
# positive m 
Weyl4_l4_m1_data = []
Weyl4_l4_m2_data = []
Weyl4_l4_m3_data = []
Weyl4_l4_m4_data = []
# negative m  
Weyl4_l4_m1n_data = []
Weyl4_l4_m2n_data = []
Weyl4_l4_m3n_data = []
Weyl4_l4_m4n_data = []

Sr_data = []
trianglearea = []
timedata = []
Cycletimedata = []

# Definitions for Quadrature sceme 

N = 131

coord = np.loadtxt('PointDistFiles/lebedev/lebedev_%03d.txt'% N)
theta = coord[:,1]*pi/180; phi = coord[:,0]*pi/180 + pi
w = coord[:,2]

# Spinweighted spherical harmonics s = -2 
# =========================================
# L = 2
# =========================================
# m = 0 zero 
sY_l2_m0 = np.sqrt(15./(32.*np.pi))*np.sin(theta)**2
# positive m values : 1,2
sY_l2_m1 = np.exp(1j*phi)*np.sqrt(5./(np.pi))*np.cos(theta/2)**3*np.sin(theta/2)
sY_l2_m2 = 1./2.*np.exp(2*1j*phi)*np.sqrt(5./(np.pi))*np.cos(theta/2)**4
# negative m values :-1,-2
sY_l2_m1n = 1./2.*np.exp(-1j*phi)*np.sqrt(5./(np.pi))*np.sin(theta/2)**2*np.sin(theta)
sY_l2_m2n = 1./2.*np.exp(-2*1j*phi)*np.sqrt(5./(np.pi))*np.sin(theta/2)**4
# =========================================
# L = 3 
# =========================================
# m = 0 zero 
sY_l3_m0 = 1./4.*np.sqrt(105./(2.*np.pi))*np.cos(theta)*np.sin(theta)**2
# positive m values : 1, 2 ... 
sY_l3_m1 = 1./2.*np.exp(1j*phi)*np.sqrt(35./(2.*np.pi))*np.cos(theta/2)**3*(-1+3*np.cos(theta))*np.sin(theta/2)
sY_l3_m2 = 1./2.*np.exp(2*1j*phi)*np.sqrt(7./(np.pi))*np.cos(theta/2)**4*(-2+3*np.cos(theta))	
sY_l3_m3 = -np.exp(3*1j*phi)*np.sqrt(21./(2*np.pi))*np.cos(theta/2)**5*np.sin(theta/2)
# negative m values : -1, -2 ... 
sY_l3_m1n = 1./2.*np.exp(-1j*phi)*np.sqrt(35./(2*np.pi))*np.cos(theta/2)*(1+3*np.cos(theta))*np.sin(theta/2)**3
sY_l3_m2n = 1./2.*np.exp(-2*1j*phi)*np.sqrt(7./(np.pi))*(2+3*np.cos(theta))*np.sin(theta/2)**4
sY_l3_m3n = 1./2.*np.exp(-3*1j*phi)*np.sqrt(21./(2*np.pi))*np.sin(theta/2)**4*np.sin(theta)
# =========================================
# L = 4
# =========================================
sY_l4_m0 = 3./16.*np.sqrt(5./(2.*np.pi))*(5.+7.*np.cos(2.*theta))*np.sin(theta)**2
# positive m = 1,2,3,4
sY_l4_m1 = 3/(2*np.sqrt(2*np.pi))*np.exp(1j*phi)*np.cos(theta/2)**3*(6-7*np.cos(theta)+7*np.cos(2*theta))*np.sin(theta/2)
sY_l4_m2 = 3/(4*np.sqrt(np.pi))*np.exp(2.*1j*phi)*(9-14*np.cos(theta)+7*np.cos(2*theta))*np.cos(theta/2)**4
sY_l4_m3 = -3.*np.sqrt(7./(2.*np.pi))*np.cos(theta/2)**5*np.exp(3.*1j*phi)*(-1+2*np.cos(theta))*np.sin(theta/2)
sY_l4_m4 = 3.*np.exp(4*1j*phi)*np.sqrt(7/np.pi)*np.cos(theta/2)**6*np.sin(theta/2)**2
# negative m = -1,-2,-3,-4
sY_l4_m1n = 3./(2.*np.sqrt(2*np.pi))*np.exp(-1j*phi)*np.cos(theta/2)*(6+7*np.cos(theta)+7*np.cos(2*theta))*np.sin(theta/2)**3
sY_l4_m2n = 3./(4.*np.sqrt(np.pi))*np.exp(-2*1j*phi)*(9+14*np.cos(theta)+7*np.cos(2*theta))*np.sin(theta/2)**4
sY_l4_m3n = 3.*np.sqrt(7./(2.*np.pi))*np.exp(-3*1j*phi)*(1+2*np.cos(theta))*np.sin(theta/2)**5*np.cos(theta/2)
sY_l4_m4n = 3./4.*np.exp(-4*1j*phi)*np.sqrt(7/np.pi)*np.sin(theta/2)**4*np.sin(theta)**2
	
# ==============================
#         Loop over all frames 
# ==============================
for i in ds:

	startCycle = time.time()

	i.print_stats() 
	center =  (i.domain_right_edge/2)

	start = time.time()

# ==================================================
	# Initalising  

	print ("Integrating ... ")	
	
	Sr = 0 + 1j*0	
	# l = 2  
	Weyl4_l2_m0 = 0 + 1j*0
	# positive m  
	Weyl4_l2_m1 = 0 + 1j*0
	Weyl4_l2_m2 = 0 + 1j*0
	# negative m  
	Weyl4_l2_m1n = 0 + 1j*0 
	Weyl4_l2_m2n = 0 + 1j*0

	# l = 3  
	Weyl4_l3_m0 = 0 + 1j*0
	# positive m  
	Weyl4_l3_m1 = 0 + 1j*0
	Weyl4_l3_m2 = 0 + 1j*0
	Weyl4_l3_m3 = 0 + 1j*0

	# negative m  
	Weyl4_l3_m1n = 0 + 1j*0
	Weyl4_l3_m2n = 0 + 1j*0
	Weyl4_l3_m3n = 0 + 1j*0

	# l = 4  
	Weyl4_l4_m0 = 0 + 1j*0
	# positive m  
	Weyl4_l4_m1 = 0 + 1j*0
	Weyl4_l4_m2 = 0 + 1j*0
	Weyl4_l4_m3 = 0 + 1j*0
	Weyl4_l4_m4 = 0 + 1j*0
	# negative m  
	Weyl4_l4_m1n = 0 + 1j*0
	Weyl4_l4_m2n = 0 + 1j*0
	Weyl4_l4_m3n = 0 + 1j*0
	Weyl4_l4_m4n = 0 + 1j*0

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
		

		Weyl4_l2_m0 += 4*pi*w[k]*np.conjugate(sY_l2_m0[k])*Weyl4*Rad
		# positive m  
		Weyl4_l2_m1 += 4*pi*w[k]*np.conjugate(sY_l2_m1[k])*Weyl4*Rad
		Weyl4_l2_m2 += 4*pi*w[k]*np.conjugate(sY_l2_m2[k])*Weyl4*Rad
		# negative m  
		Weyl4_l2_m1n += 4*pi*w[k]*np.conjugate(sY_l2_m1n[k])*Weyl4*Rad
		Weyl4_l2_m2n += 4*pi*w[k]*np.conjugate(sY_l2_m2n[k])*Weyl4*Rad
		# l = 3  
		Weyl4_l3_m0 += 4*pi*w[k]*np.conjugate(sY_l3_m0[k])*Weyl4*Rad
		# positive m  
		Weyl4_l3_m1 += 4*pi*w[k]*np.conjugate(sY_l3_m1[k])*Weyl4*Rad
		Weyl4_l3_m2 += 4*pi*w[k]*np.conjugate(sY_l3_m2[k])*Weyl4*Rad
		Weyl4_l3_m3 += 4*pi*w[k]*np.conjugate(sY_l3_m3[k])*Weyl4*Rad
		# negative m  
		Weyl4_l3_m1n += 4*pi*w[k]*np.conjugate(sY_l3_m1n[k])*Weyl4*Rad
		Weyl4_l3_m2n += 4*pi*w[k]*np.conjugate(sY_l3_m2n[k])*Weyl4*Rad
		Weyl4_l3_m3n += 4*pi*w[k]*np.conjugate(sY_l3_m3n[k])*Weyl4*Rad
		# l = 4  
		Weyl4_l4_m0 += 4*pi*w[k]*np.conjugate(sY_l4_m0[k])*Weyl4*Rad
		# positive m  
		Weyl4_l4_m1 += 4*pi*w[k]*np.conjugate(sY_l4_m1[k])*Weyl4*Rad
		Weyl4_l4_m2 += 4*pi*w[k]*np.conjugate(sY_l4_m2[k])*Weyl4*Rad
		Weyl4_l4_m3 += 4*pi*w[k]*np.conjugate(sY_l4_m3[k])*Weyl4*Rad
		Weyl4_l4_m4 += 4*pi*w[k]*np.conjugate(sY_l4_m4[k])*Weyl4*Rad
		# negative m  
		Weyl4_l4_m1n += 4*pi*w[k]*np.conjugate(sY_l4_m1n[k])*Weyl4*Rad
		Weyl4_l4_m2n += 4*pi*w[k]*np.conjugate(sY_l4_m2n[k])*Weyl4*Rad
		Weyl4_l4_m3n += 4*pi*w[k]*np.conjugate(sY_l4_m3n[k])*Weyl4*Rad
		Weyl4_l4_m4n += 4*pi*w[k]*np.conjugate(sY_l4_m4n[k])*Weyl4*Rad

# ==================================================
# DATA WRITEOUT
# ==================================================


	
	# l = 2  
	Weyl4_l2_m0_data.append(Weyl4_l2_m0)
	# positive m  
	Weyl4_l2_m1_data.append(Weyl4_l2_m1)  
	Weyl4_l2_m2_data.append(Weyl4_l2_m2)  
	# negative m  
	Weyl4_l2_m1n_data.append(Weyl4_l2_m1n)  
	Weyl4_l2_m2n_data.append(Weyl4_l2_m2n)  


	# l = 3  
	Weyl4_l3_m0_data.append(Weyl4_l3_m0)
	# positive m  
	Weyl4_l3_m1_data.append(Weyl4_l3_m1)  
	Weyl4_l3_m2_data.append(Weyl4_l3_m2)  
	Weyl4_l3_m3_data.append(Weyl4_l3_m3) 
	# negative m  
	Weyl4_l3_m1n_data.append(Weyl4_l3_m1n)  
	Weyl4_l3_m2n_data.append(Weyl4_l3_m2n)  
	Weyl4_l3_m3n_data.append(Weyl4_l3_m3n)  

	# l = 4  
	Weyl4_l4_m0_data.append(Weyl4_l4_m0)
	# positive m  
	Weyl4_l4_m1_data.append(Weyl4_l4_m1)  
	Weyl4_l4_m2_data.append(Weyl4_l4_m2)  
	Weyl4_l4_m3_data.append(Weyl4_l4_m3) 
	Weyl4_l4_m4_data.append(Weyl4_l4_m4) 
	# negative m  
	Weyl4_l4_m1n_data.append(Weyl4_l4_m1n) 
	Weyl4_l4_m2n_data.append(Weyl4_l4_m2n)  
	Weyl4_l4_m3n_data.append(Weyl4_l4_m3n) 
	Weyl4_l4_m4n_data.append(Weyl4_l4_m4n)

	# time
	timedata.append(i.current_time)

	np.savetxt('time.out',timedata)
	# l = 2 
	np.savetxt('Weyl4_l2_m0_data.out',Weyl4_l2_m0_data)
	np.savetxt('Weyl4_l2_m1_data.out',Weyl4_l2_m1_data)
	np.savetxt('Weyl4_l2_m2_data.out',Weyl4_l2_m2_data)
	np.savetxt('Weyl4_l2_m1n_data.out',Weyl4_l2_m1n_data)
	np.savetxt('Weyl4_l2_m2n_data.out',Weyl4_l2_m2n_data)
	# l = 3 
	np.savetxt('Weyl4_l3_m0_data.out',Weyl4_l3_m0_data)
	np.savetxt('Weyl4_l3_m1_data.out',Weyl4_l3_m1_data)
	np.savetxt('Weyl4_l3_m2_data.out',Weyl4_l3_m2_data)
	np.savetxt('Weyl4_l3_m3_data.out',Weyl4_l3_m3_data)
	np.savetxt('Weyl4_l3_m1n_data.out',Weyl4_l3_m1n_data)
	np.savetxt('Weyl4_l3_m2n_data.out',Weyl4_l3_m2n_data)
	np.savetxt('Weyl4_l3_m3n_data.out',Weyl4_l3_m3n_data)	
	# l = 4
	np.savetxt('Weyl4_l4_m0_data.out',Weyl4_l4_m0_data)
	np.savetxt('Weyl4_l4_m1_data.out',Weyl4_l4_m1_data)
	np.savetxt('Weyl4_l4_m2_data.out',Weyl4_l4_m2_data)
	np.savetxt('Weyl4_l4_m3_data.out',Weyl4_l4_m3_data)
	np.savetxt('Weyl4_l4_m4_data.out',Weyl4_l4_m4_data)
	np.savetxt('Weyl4_l4_m1n_data.out',Weyl4_l4_m1n_data)
	np.savetxt('Weyl4_l4_m2n_data.out',Weyl4_l4_m2n_data)
	np.savetxt('Weyl4_l4_m3n_data.out',Weyl4_l4_m3n_data)
	np.savetxt('Weyl4_l4_m4n_data.out',Weyl4_l4_m4n_data)

	mid1 = time.time()
	print("time to Integrate= ", mid1-start)

# ==================================================
# Performance 
# ==================================================
	Cycletime = abs(time.time()-startCycle)
	Cycletimedata.append(Cycletime)
	np.savetxt('Cycletime.out',Cycletimedata)
# ==================================================



print("FINISHED !!! ")


