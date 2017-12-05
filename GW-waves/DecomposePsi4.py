import yt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
import time
from numpy import sin,cos,sqrt 
 
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
loc = '../../GRChombo/GRChombo_1.0/S1x0000*'

def csc(x):
	return 1./cos(x)

def sec(x):
	return 1./sin(x)

def Sqrt(x):
	return np.sqrt(x)
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
# l = 5 
Weyl4_l5_m0_data = []
# positive m 
Weyl4_l5_m1_data = []
Weyl4_l5_m2_data = []
Weyl4_l5_m3_data = []
Weyl4_l5_m4_data = []
Weyl4_l5_m5_data = []
# negative m  
Weyl4_l5_m1n_data = []
Weyl4_l5_m2n_data = []
Weyl4_l5_m3n_data = []
Weyl4_l5_m4n_data = []
Weyl4_l5_m5n_data = []
# l = 6 
Weyl4_l6_m0_data = []
# positive m 
Weyl4_l6_m1_data = []
Weyl4_l6_m2_data = []
Weyl4_l6_m3_data = []
Weyl4_l6_m4_data = []
Weyl4_l6_m5_data = []
Weyl4_l6_m6_data = []
# negative m  
Weyl4_l6_m1n_data = []
Weyl4_l6_m2n_data = []
Weyl4_l6_m3n_data = []
Weyl4_l6_m4n_data = []
Weyl4_l6_m5n_data = []
Weyl4_l6_m6n_data = []

trianglearea = []
timedata = []
Cycletimedata = []

Euler = np.e 
Pi = np.pi
I = 1j

print (csc(1)-1/cos(1))

@derived_field(name = "rr", units = "cm")
def _rr(field, data):
    center = i.domain_right_edge/2
    x = (data["x"]-center[0])
    y = (data["y"]-center[1])
    z = (data["z"]-center[2])
    return np.sqrt(x*x+y*y+z*z)

# ==============================
#         Loop over all frames 
# ==============================
for i in ds:

	startCycle = time.time()

	i.print_stats() 
	center =  (i.domain_right_edge/2)
	surf = i.surface(i.all_data(),"rr",Rad)
	start = time.time()
	trianglearea = []
	triangles = surf.triangles
	totalarea = 0

# ==================================================
# Calulcating area of every triangle on the surface

	for j in range(len(triangles)):		
		a = triangles[j,0] - triangles[j,1]
		b = triangles[j,0] - triangles[j,2]	
		triarea = 1./2*np.linalg.norm(np.cross(a,b))
		trianglearea.append(triarea)
		totalarea += triarea
		
	print ("area of sphere = ", totalarea)
	mid1 = time.time()
	print("time to calculate triangles = ", mid1-start)

# ==================================================
# Getting position data

	x = (surf["x"]-float(center[0]))
	y = (surf["y"]-float(center[1]))
	z = (surf["z"]-float(center[2]))
	rrho = np.sqrt(x*x+y*y)
	rr = np.sqrt(x*x+y*y+z*z)
	sinth = rrho/rr
	costh = z/rr	
	 	
	theta = np.arccos(costh)
	phi = np.where(y > 0,np.arccos((x/rrho)),-np.arccos((x/rrho)))
	
	print(np.linalg.norm(np.sin(theta)-sinth))

	start = time.time()	

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
	sY_l4_m0 = (3*sqrt(5/(2.*Pi))*(5 + 7*cos(2*theta))*sin(theta)**2)/16. 
	# positive m = 1,2,3,4
	sY_l4_m1 = (3*Euler**(I*phi)*cos(theta/2.)**3*(6 - 7*cos(theta) + 7*cos(2*theta))*sin(theta/2.))/(2.*Sqrt(2*Pi)) 
	sY_l4_m2 = (3*Euler**(2*I*phi)*cos(theta/2.)**4*(9 - 14*cos(theta) + 7*cos(2*theta)))/(4.*Sqrt(Pi))
	sY_l4_m3 = -3*Euler**(3*I*phi)*Sqrt(7/(2.*Pi))*cos(theta/2.)**5*(-1 + 2*cos(theta))*sin(theta/2.)
	sY_l4_m4 = (3*Euler**(4*I*phi)*Sqrt(7/Pi)*csc(theta/2.)**4*sin(theta)**6)/64.
 	# negative m = -1,-2,-3,-4
	sY_l4_m1n = (3*cos(theta/2.)*(6 + 7*cos(theta) + 7*cos(2*theta))*sin(theta/2.)**3)/(2.*Euler**(I*phi)*Sqrt(2*Pi))
	sY_l4_m2n = (3*(9 + 14*cos(theta) + 7*cos(2*theta))*sin(theta/2.)**4)/(4.*Euler**(2*I*phi)*Sqrt(Pi))
	sY_l4_m3n = (3*Sqrt(7/(2.*Pi))*cos(theta/2.)*(1 + 2*cos(theta))*sin(theta/2.)**5)/Euler**(3*I*phi)
	sY_l4_m4n = (3*Sqrt(7/Pi)*sin(theta/2.)**4*sin(theta)**2)/(4.*Euler**(4*I*phi))
	
# =========================================
# L = 5  
# =========================================
	sY_l5_m0 = (Sqrt(1155/(2.*Pi))*(5*cos(theta) + 3*cos(3*theta))*sin(theta)**2)/32. 
	# positive m = 1,2,3,4,5
	sY_l5_m1 = (Euler**(I*phi)*Sqrt(77/Pi)*cos(theta/2.)**3*(-14 + 33*cos(theta) - 18*cos(2*theta) + 15*cos(3*theta))*sin(theta/2.))/16.
	sY_l5_m2 = (Euler**(2*I*phi)*Sqrt(11/Pi)*cos(theta/2.)**4*(-32 + 57*cos(theta) - 36*cos(2*theta) + 15*cos(3*theta)))/8.
	sY_l5_m3 = -(Euler**(3*I*phi)*Sqrt(33/(2.*Pi))*cos(theta/2.)**5*(17 - 24*cos(theta) + 15*cos(2*theta))*sin(theta/2.))/4.
	sY_l5_m4 = Euler**(4*I*phi)*Sqrt(33/Pi)*cos(theta/2.)**6*(-2 + 5*cos(theta))*sin(theta/2.)**2
	sY_l5_m5 = -(Euler**(5*I*phi)*Sqrt(330/Pi)*cos(theta/2.)**7*sin(theta/2.)**3)
	# negative m = -1,-2,-3,-4,-5
	sY_l5_m1n = (Sqrt(77/Pi)*cos(theta/2.)*(14 + 33*cos(theta) + 18*cos(2*theta) + 15*cos(3*theta))*sin(theta/2.)**3)/(16.*Euler**(I*phi))
	sY_l5_m2n = (Sqrt(11/Pi)*(32 + 57*cos(theta) + 36*cos(2*theta) + 15*cos(3*theta))*sin(theta/2.)**4)/(8.*Euler**(2*I*phi))
	sY_l5_m3n = (Sqrt(33/(2.*Pi))*cos(theta/2.)*(17 + 24*cos(theta) + 15*cos(2*theta))*sin(theta/2.)**5)/(4.*Euler**(3*I*phi))
	sY_l5_m4n = (Sqrt(33/Pi)*cos(theta/2.)**2*(2 + 5*cos(theta))*sin(theta/2.)**6)/Euler**(4*I*phi)
	sY_l5_m5n = (Sqrt(330/Pi)*cos(theta/2.)**3*sin(theta/2.)**7)/Euler**(5*I*phi) 
	
# =========================================
# L = 6  
# =========================================
	sY_l6_m0 = (Sqrt(1365/Pi)*(35 + 60*cos(2*theta) + 33*cos(4*theta))*sin(theta)**2)/512. 
	# positive m = 1,2,3,4,5,6
	sY_l6_m1 = (Euler**(I*phi)*Sqrt(65/(2.*Pi))*cos(theta/2.)**3*(161 - 252*cos(theta) + 252*cos(2*theta) - 132*cos(3*theta) + 99*cos(4*theta))*sin(theta/2.))/64.
	sY_l6_m2 = (Euler**(2*I*phi)*Sqrt(13/Pi)*cos(theta/2.)**4*(1709 - 3096*cos(theta) + 2340*cos(2*theta) - 1320*cos(3*theta) + 495*cos(4*theta)))/256.
	sY_l6_m3 = (-3*Euler**(3*I*phi)*Sqrt(13/Pi)*cos(theta/2.)**5*(-98 + 185*cos(theta) - 110*cos(2*theta) + 55*cos(3*theta))*sin(theta/2.))/32.
	sY_l6_m4 = (Euler**(4*I*phi)*Sqrt(195/(2.*Pi))*cos(theta/2.)**6*(35 - 44*cos(theta) + 33*cos(2*theta))*sin(theta/2.)**2)/8.
	sY_l6_m5 = -(Euler**(5*I*phi)*Sqrt(2145/Pi)*cos(theta/2.)**7*(-1 + 3*cos(theta))*sin(theta/2.)**3)/2.
	sY_l6_m6 = (3*Euler**(6*I*phi)*Sqrt(715/Pi)*csc(theta/2.)**4*sin(theta)**8)/512.
 	# negative m = -1,-2,-3,-4,-5,-6
	sY_l6_m1n = (Sqrt(65/(2.*Pi))*cos(theta/2.)*(161 + 252*cos(theta) + 252*cos(2*theta) + 132*cos(3*theta) + 99*cos(4*theta))*sin(theta/2.)**3)/(64.*Euler**(I*phi)) 
	sY_l6_m2n = (Sqrt(13/Pi)*(1709 + 3096*cos(theta) + 2340*cos(2*theta) + 1320*cos(3*theta) + 495*cos(4*theta))*sin(theta/2.)**4)/(256.*Euler**(2*I*phi))
	sY_l6_m3n = (3*Sqrt(13/Pi)*cos(theta/2.)*(98 + 185*cos(theta) + 110*cos(2*theta) + 55*cos(3*theta))*sin(theta/2.)**5)/(32.*Euler**(3*I*phi))
	sY_l6_m4n = (Sqrt(195/(2.*Pi))*cos(theta/2.)**2*(35 + 44*cos(theta) + 33*cos(2*theta))*sin(theta/2.)**6)/(8.*Euler**(4*I*phi))
	sY_l6_m5n = (Sqrt(2145/Pi)*cos(theta/2.)**3*(1 + 3*cos(theta))*sin(theta/2.)**7)/(2.*Euler**(5*I*phi))
	sY_l6_m6n = (3*Sqrt(715/Pi)*sin(theta/2.)**4*sin(theta)**4)/(32.*Euler**(6*I*phi))
	
	print("time to define Spher Harms= ", time.time()-start)

# ==================================================
	# Getting Weyl4 
	# This is the computiational most intensive task 
	
	start = time.time()	

	Weyl4 = surf["Weyl1"] + 1j*surf["Weyl2"]
	
	print("time to extract Weyl_4 = ", time.time()-start)

# ==================================================
	# Initalising  

	print ("Integrating ... ")	
	
	# l = 2  
	Weyl4_l2_m0 = 0
	# positive m  
	Weyl4_l2_m1 = 0 
	Weyl4_l2_m2 = 0
	# negative m  
	Weyl4_l2_m1n = 0 
	Weyl4_l2_m2n = 0

	# l = 3  
	Weyl4_l3_m0 = 0
	# positive m  
	Weyl4_l3_m1 = 0 
	Weyl4_l3_m2 = 0 
	Weyl4_l3_m3 = 0
	# negative m  
	Weyl4_l3_m1n = 0 
	Weyl4_l3_m2n = 0 
	Weyl4_l3_m3n = 0 

	# l = 4  
	Weyl4_l4_m0 = 0
	# positive m  
	Weyl4_l4_m1 = 0 
	Weyl4_l4_m2 = 0 
	Weyl4_l4_m3 = 0
	Weyl4_l4_m4 = 0
	# negative m  
	Weyl4_l4_m1n = 0 
	Weyl4_l4_m2n = 0 
	Weyl4_l4_m3n = 0 
	Weyl4_l4_m4n = 0

	# l = 5 
	Weyl4_l5_m0 = 0
	# positive m  
	Weyl4_l5_m1 = 0 
	Weyl4_l5_m2 = 0 
	Weyl4_l5_m3 = 0
	Weyl4_l5_m4 = 0
	Weyl4_l5_m5 = 0
	# negative m  
	Weyl4_l5_m1n = 0 
	Weyl4_l5_m2n = 0 
	Weyl4_l5_m3n = 0 
	Weyl4_l5_m4n = 0
	Weyl4_l5_m5n = 0

	# l = 6 
	Weyl4_l6_m0 = 0
	# positive m  
	Weyl4_l6_m1 = 0 
	Weyl4_l6_m2 = 0 
	Weyl4_l6_m3 = 0
	Weyl4_l6_m4 = 0
	Weyl4_l6_m5 = 0
	Weyl4_l6_m6 = 0
	# negative m  
	Weyl4_l6_m1n = 0 
	Weyl4_l6_m2n = 0 
	Weyl4_l6_m3n = 0 
	Weyl4_l6_m4n = 0
	Weyl4_l6_m5n = 0
	Weyl4_l6_m6n = 0



	start = time.time()
	
	index = 0 
	for k in range(len(trianglearea)):
		Weyl4_l2_m0 += trianglearea[k]*np.conjugate(sY_l2_m0[k])*Weyl4[k]/Rad
		# positive m  
		Weyl4_l2_m1 += trianglearea[k]*np.conjugate(sY_l2_m1[k])*Weyl4[k]/Rad
		Weyl4_l2_m2 += trianglearea[k]*np.conjugate(sY_l2_m2[k])*Weyl4[k]/Rad
		# negative m  
		Weyl4_l2_m1n += trianglearea[k]*np.conjugate(sY_l2_m1n[k])*Weyl4[k]/Rad
		Weyl4_l2_m2n += trianglearea[k]*np.conjugate(sY_l2_m2n[k])*Weyl4[k]/Rad
		# l = 3  
		Weyl4_l3_m0 += trianglearea[k]*np.conjugate(sY_l3_m0[k])*Weyl4[k]/Rad
		# positive m  
		Weyl4_l3_m1 += trianglearea[k]*np.conjugate(sY_l3_m1[k])*Weyl4[k]/Rad
		Weyl4_l3_m2 += trianglearea[k]*np.conjugate(sY_l3_m2[k])*Weyl4[k]/Rad
		Weyl4_l3_m3 += trianglearea[k]*np.conjugate(sY_l3_m3[k])*Weyl4[k]/Rad
		# negative m  
		Weyl4_l3_m1n += trianglearea[k]*np.conjugate(sY_l3_m1n[k])*Weyl4[k]/Rad
		Weyl4_l3_m2n += trianglearea[k]*np.conjugate(sY_l3_m2n[k])*Weyl4[k]/Rad
		Weyl4_l3_m3n += trianglearea[k]*np.conjugate(sY_l3_m3n[k])*Weyl4[k]/Rad
		# l = 4  
		Weyl4_l4_m0 += trianglearea[k]*np.conjugate(sY_l4_m0[k])*Weyl4[k]/Rad
		# positive m  
		Weyl4_l4_m1 += trianglearea[k]*np.conjugate(sY_l4_m1[k])*Weyl4[k]/Rad
		Weyl4_l4_m2 += trianglearea[k]*np.conjugate(sY_l4_m2[k])*Weyl4[k]/Rad
		Weyl4_l4_m3 += trianglearea[k]*np.conjugate(sY_l4_m3[k])*Weyl4[k]/Rad
		Weyl4_l4_m4 += trianglearea[k]*np.conjugate(sY_l4_m4[k])*Weyl4[k]/Rad
		# negative m  
		Weyl4_l4_m1n += trianglearea[k]*np.conjugate(sY_l4_m1n[k])*Weyl4[k]/Rad
		Weyl4_l4_m2n += trianglearea[k]*np.conjugate(sY_l4_m2n[k])*Weyl4[k]/Rad
		Weyl4_l4_m3n += trianglearea[k]*np.conjugate(sY_l4_m3n[k])*Weyl4[k]/Rad
		Weyl4_l4_m4n += trianglearea[k]*np.conjugate(sY_l4_m4n[k])*Weyl4[k]/Rad
		# l = 5  
		Weyl4_l5_m0 += trianglearea[k]*np.conjugate(sY_l5_m0[k])*Weyl4[k]/Rad
		# positive m  
		Weyl4_l5_m1 += trianglearea[k]*np.conjugate(sY_l5_m1[k])*Weyl4[k]/Rad
		Weyl4_l5_m2 += trianglearea[k]*np.conjugate(sY_l5_m2[k])*Weyl4[k]/Rad
		Weyl4_l5_m3 += trianglearea[k]*np.conjugate(sY_l5_m3[k])*Weyl4[k]/Rad
		Weyl4_l5_m4 += trianglearea[k]*np.conjugate(sY_l5_m4[k])*Weyl4[k]/Rad
		Weyl4_l5_m5 += trianglearea[k]*np.conjugate(sY_l5_m5[k])*Weyl4[k]/Rad
		# negative m  
		Weyl4_l5_m1n += trianglearea[k]*np.conjugate(sY_l5_m1n[k])*Weyl4[k]/Rad
		Weyl4_l5_m2n += trianglearea[k]*np.conjugate(sY_l5_m2n[k])*Weyl4[k]/Rad
		Weyl4_l5_m3n += trianglearea[k]*np.conjugate(sY_l5_m3n[k])*Weyl4[k]/Rad
		Weyl4_l5_m4n += trianglearea[k]*np.conjugate(sY_l5_m4n[k])*Weyl4[k]/Rad
		Weyl4_l5_m5n += trianglearea[k]*np.conjugate(sY_l5_m5n[k])*Weyl4[k]/Rad
		# l = 6
		Weyl4_l6_m0 += trianglearea[k]*np.conjugate(sY_l6_m0[k])*Weyl4[k]/Rad
		# positive m  
		Weyl4_l6_m1 += trianglearea[k]*np.conjugate(sY_l6_m1[k])*Weyl4[k]/Rad
		Weyl4_l6_m2 += trianglearea[k]*np.conjugate(sY_l6_m2[k])*Weyl4[k]/Rad
		Weyl4_l6_m3 += trianglearea[k]*np.conjugate(sY_l6_m3[k])*Weyl4[k]/Rad
		Weyl4_l6_m4 += trianglearea[k]*np.conjugate(sY_l6_m4[k])*Weyl4[k]/Rad
		Weyl4_l6_m5 += trianglearea[k]*np.conjugate(sY_l6_m5[k])*Weyl4[k]/Rad
		Weyl4_l6_m6 += trianglearea[k]*np.conjugate(sY_l6_m6[k])*Weyl4[k]/Rad
		# negative m  
		Weyl4_l6_m1n += trianglearea[k]*np.conjugate(sY_l6_m1n[k])*Weyl4[k]/Rad
		Weyl4_l6_m2n += trianglearea[k]*np.conjugate(sY_l6_m2n[k])*Weyl4[k]/Rad
		Weyl4_l6_m3n += trianglearea[k]*np.conjugate(sY_l6_m3n[k])*Weyl4[k]/Rad
		Weyl4_l6_m4n += trianglearea[k]*np.conjugate(sY_l6_m4n[k])*Weyl4[k]/Rad
		Weyl4_l6_m5n += trianglearea[k]*np.conjugate(sY_l6_m5n[k])*Weyl4[k]/Rad
		Weyl4_l6_m6n += trianglearea[k]*np.conjugate(sY_l6_m6n[k])*Weyl4[k]/Rad

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

	# l = 5  
	Weyl4_l5_m0_data.append(Weyl4_l5_m0)
	# positive m  
	Weyl4_l5_m1_data.append(Weyl4_l5_m1)  
	Weyl4_l5_m2_data.append(Weyl4_l5_m2)  
	Weyl4_l5_m3_data.append(Weyl4_l5_m3) 
	Weyl4_l5_m4_data.append(Weyl4_l5_m4) 
	Weyl4_l5_m5_data.append(Weyl4_l5_m5) 
	# negative m  
	Weyl4_l5_m1n_data.append(Weyl4_l5_m1n) 
	Weyl4_l5_m2n_data.append(Weyl4_l5_m2n)  
	Weyl4_l5_m3n_data.append(Weyl4_l5_m3n) 
	Weyl4_l5_m4n_data.append(Weyl4_l5_m4n)
	Weyl4_l5_m5n_data.append(Weyl4_l5_m5n)

	# l = 6 
	Weyl4_l6_m0_data.append(Weyl4_l6_m0)
	# positive m  
	Weyl4_l6_m1_data.append(Weyl4_l6_m1)  
	Weyl4_l6_m2_data.append(Weyl4_l6_m2)  
	Weyl4_l6_m3_data.append(Weyl4_l6_m3) 
	Weyl4_l6_m4_data.append(Weyl4_l6_m4) 
	Weyl4_l6_m5_data.append(Weyl4_l6_m5) 
	Weyl4_l6_m6_data.append(Weyl4_l6_m6) 
	# negative m  
	Weyl4_l6_m1n_data.append(Weyl4_l6_m1n) 
	Weyl4_l6_m2n_data.append(Weyl4_l6_m2n)  
	Weyl4_l6_m3n_data.append(Weyl4_l6_m3n) 
	Weyl4_l6_m4n_data.append(Weyl4_l6_m4n)
	Weyl4_l6_m5n_data.append(Weyl4_l6_m5n)
	Weyl4_l6_m6n_data.append(Weyl4_l6_m6n)



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
	# l = 5
	np.savetxt('Weyl4_l5_m0_data.out',Weyl4_l5_m0_data)
	np.savetxt('Weyl4_l5_m1_data.out',Weyl4_l5_m1_data)
	np.savetxt('Weyl4_l5_m2_data.out',Weyl4_l5_m2_data)
	np.savetxt('Weyl4_l5_m3_data.out',Weyl4_l5_m3_data)
	np.savetxt('Weyl4_l5_m4_data.out',Weyl4_l5_m4_data)
	np.savetxt('Weyl4_l5_m5_data.out',Weyl4_l5_m5_data)
	np.savetxt('Weyl4_l5_m1n_data.out',Weyl4_l5_m1n_data)
	np.savetxt('Weyl4_l5_m2n_data.out',Weyl4_l5_m2n_data)
	np.savetxt('Weyl4_l5_m3n_data.out',Weyl4_l5_m3n_data)
	np.savetxt('Weyl4_l5_m4n_data.out',Weyl4_l5_m4n_data)
	np.savetxt('Weyl4_l5_m5n_data.out',Weyl4_l5_m5n_data)
	# l = 6
	np.savetxt('Weyl4_l6_m0_data.out',Weyl4_l6_m0_data)
	np.savetxt('Weyl4_l6_m1_data.out',Weyl4_l6_m1_data)
	np.savetxt('Weyl4_l6_m2_data.out',Weyl4_l6_m2_data)
	np.savetxt('Weyl4_l6_m3_data.out',Weyl4_l6_m3_data)
	np.savetxt('Weyl4_l6_m4_data.out',Weyl4_l6_m4_data)
	np.savetxt('Weyl4_l6_m5_data.out',Weyl4_l6_m5_data)
	np.savetxt('Weyl4_l6_m6_data.out',Weyl4_l6_m6_data)
	np.savetxt('Weyl4_l6_m1n_data.out',Weyl4_l6_m1n_data)
	np.savetxt('Weyl4_l6_m2n_data.out',Weyl4_l6_m2n_data)
	np.savetxt('Weyl4_l6_m3n_data.out',Weyl4_l6_m3n_data)
	np.savetxt('Weyl4_l6_m4n_data.out',Weyl4_l6_m4n_data)
	np.savetxt('Weyl4_l6_m5n_data.out',Weyl4_l6_m5n_data)
	np.savetxt('Weyl4_l6_m6n_data.out',Weyl4_l6_m6n_data)

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


