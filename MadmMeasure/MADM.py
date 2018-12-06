import yt
from yt import derived_field
import numpy as np
from numpy.linalg import inv
import time
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import pi
yt.enable_parallelism()

# ==============================================================================
#			Spherical Horizon finder
#
#  Assumes (Obviously) spherical symmetry !! Gives bad results if not !
#  Input = Black hole center
#
#  For valid calculation of Mass one has to slice along x-Axis !!
#       Output :-Black Hole mass (with Error estimate)
#               -Black Hole radius (with Error estimate)
# ==============================================================================

#Loading dataset
ds = yt.load('../../KerrBH_*.hdf5')

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])

Rad = 250.

NInt = 27
coord = np.loadtxt('PointDistFiles/lebedev/lebedev_%03d.txt'% NInt)
theta = coord[:,1]*pi/180; phi = coord[:,0]*pi/180 + pi
w = coord[:,2]

NumInt = len(theta)
print("We are using", NumInt, " point over the sphere")

# =================================
#  INPUT #
# =================================
BHcenter = [center,center,center]
# ================================
# ================================

time_data = []
CycleData = []
MADM_data = []
J1ADM_data = []
J2ADM_data = []
J3ADM_data = []

# used interpolation quality
Quality = 'quadratic'

for i in ds:

	start = time.time()


	# Gradients calculated by YT

	i.add_gradient_fields(('chombo', u'chi'))
	i.add_gradient_fields(('chombo', u'h11'))
	i.add_gradient_fields(('chombo', u'h12'))
	i.add_gradient_fields(('chombo', u'h13'))
	i.add_gradient_fields(('chombo', u'h22'))
	i.add_gradient_fields(('chombo', u'h23'))
	i.add_gradient_fields(('chombo', u'h33'))


# ========================================================
#	Definitions
# ========================================================


	print ("Preparing calculation ... ")

	counter = 0

	MADM = 0

	JADM = np.zeros(3)
	PADM = np.zeros(3)

	for (k,x) in enumerate(phi):
		phi_var = phi[k]
		theta_var = theta[k]
		x1 = Rad*np.cos(phi_var)*np.sin(theta_var)+float(centerXYZ[0])
		y1 = Rad*np.sin(phi_var)*np.sin(theta_var)+float(centerXYZ[1])
		z1 = Rad*np.cos(theta_var)+float(centerXYZ[2])
		ptnEx = [x1,y1,z1]
		c = i.point(ptnEx)

		BHcenterYT = i.arr(BHcenter, "code_length")


	# Yt gives data unsorted

		x = np.array(c["x"]- BHcenterYT[0])
		y = np.array(c["y"]- BHcenterYT[1])
		z = np.array(c["z"]- BHcenterYT[2])

		r = np.sqrt(x*x+y*y+z*z)
		srt = np.argsort(r)

	# Importing Coordinates

		x = np.array(c["x"][srt]- BHcenterYT[0])
		y = np.array(c["y"][srt]- BHcenterYT[1])
		z = np.array(c["z"][srt]- BHcenterYT[2])

		pos = []

		pos.append(x)
		pos.append(y)
		pos.append(z)

		rho2 = x*x + y*y
		rho = np.sqrt(rho2)

		r2 = rho2 + z*z
		r = np.sqrt(r2)

		cosphi = x/rho
		sinphi = y/rho
		sintheta = rho/r
		costheta = z/r

		chi = np.array(c["chi"][srt]);
		chi2 = chi*chi


		g11 = np.array(c["h11"][srt]);
		g12 = np.array(c["h12"][srt]);
		g13 = np.array(c["h13"][srt]);
		g22 = np.array(c["h22"][srt]);
		g23 = np.array(c["h23"][srt]);
		g33 = np.array(c["h33"][srt]);
	
	
		g = [[],[],[]]

		g[0].append(g11)
		g[0].append(g12)
		g[0].append(g13)
		g[1].append(g12)
		g[1].append(g22)
		g[1].append(g23)
		g[2].append(g13)
		g[2].append(g23)
		g[2].append(g33)


		g = np.array(g)
		gu = np.linalg.inv(g[:,:,0])


		A11 = np.array(c["A11"][srt]);
		A12 = np.array(c["A12"][srt]);
		A13 = np.array(c["A13"][srt]);
		A22 = np.array(c["A22"][srt]);
		A23 = np.array(c["A23"][srt]);
		A33 = np.array(c["A33"][srt]);

		K = np.array(c["K"][srt]);


		A = [[],[],[]]

		A[0].append(A11)
		A[0].append(A12)
		A[0].append(A13)
		A[1].append(A12)
		A[1].append(A22)
		A[1].append(A23)
		A[2].append(A13)
		A[2].append(A23)
		A[2].append(A33)

		A = np.array(A)

# ========================================================
#	Define normal vector and invert metric
# ========================================================

		N = len(chi)
		s = np.zeros([3,N])

		for IND in range(N):
			# Definitions
			# ---------------------
			s[0,IND] = x[IND]/r[IND]
			s[1,IND] = y[IND]/r[IND]
			s[2,IND] = z[IND]/r[IND]

		

# ========================================================
#	Derivatives
# ========================================================

		dgdx = np.zeros([3,3,3,N])
		dchidx = np.zeros([3,N])

		print ("Calculating metric x Gradient ... ")
		dgdx[0,0,0,:] = c["h11_gradient_x"][srt]
		dgdx[0,1,0,:] = c["h12_gradient_x"][srt]
		dgdx[0,2,0,:] = c["h13_gradient_x"][srt]
		dgdx[1,1,0,:] = c["h22_gradient_x"][srt]
		dgdx[1,2,0,:] = c["h23_gradient_x"][srt]
		dgdx[2,2,0,:] = c["h33_gradient_x"][srt]
		dgdx[1,0,0,:] = dgdx[0,1,0,:]
		dgdx[2,0,0,:] = dgdx[0,2,0,:]
		dgdx[2,1,0,:] = dgdx[1,2,0,:]


		print ("Calculating metric y Gradient ... ")
		dgdx[0,0,1,:] = c["h11_gradient_y"][srt]
		dgdx[0,1,1,:] = c["h12_gradient_y"][srt]
		dgdx[0,2,1,:] = c["h13_gradient_y"][srt]
		dgdx[1,1,1,:] = c["h22_gradient_y"][srt]
		dgdx[1,2,1,:] = c["h23_gradient_y"][srt]
		dgdx[2,2,1,:] = c["h33_gradient_y"][srt]
		dgdx[1,0,1,:] = dgdx[0,1,1,:]
		dgdx[2,0,1,:] = dgdx[0,2,1,:]
		dgdx[2,1,1,:] = dgdx[1,2,1,:]
	

		print ("Calculating metric z Gradient ... ")
		dgdx[0,0,2,:] = c["h11_gradient_z"][srt]
		dgdx[0,1,2,:] = c["h12_gradient_z"][srt]
		dgdx[0,2,2,:] = c["h13_gradient_z"][srt]
		dgdx[1,1,2,:] = c["h22_gradient_z"][srt]
		dgdx[1,2,2,:] = c["h23_gradient_z"][srt]
		dgdx[2,2,2,:] = c["h33_gradient_z"][srt]
		dgdx[1,0,2,:] = dgdx[0,1,2,:]
		dgdx[2,0,2,:] = dgdx[0,2,2,:]
		dgdx[2,1,2,:] = dgdx[1,2,2,:]
	

		print ("Calculating chi Gradient ... ")
		dchidx[0,:] = c["chi_gradient_x"][srt]
		dchidx[1,:] = c["chi_gradient_y"][srt]
		dchidx[2,:] = c["chi_gradient_z"][srt]
	
# ========================================================
#	Calculation of Omega
#	Eq. 7.41 in shapiro book
# ========================================================

		counter += 1. 
		print ("Calculating MADM ... ")
		print ("Percentage ",counter/NumInt*100," %"  )
#		MADMInt = np.zeros(N)
		MADMInt = 0
		JADMInt = np.zeros(3)
		PADMInt = np.zeros(3)
		levi = np.zeros((3,3,3))
		Ktsr = np.zeros((3,3))
		deltaij = np.zeros((3,3))
		pstn = np.zeros(3)

		deltaij[0,0] = 1 
		deltaij[1,1] = 1 
		deltaij[2,2] = 1 

		levi[0,1,2] = 1
		levi[1,2,0] = 1
		levi[2,0,1] = 1
		levi[2,1,0] = -1
		levi[1,0,2] = -1
		levi[0,2,1] = -1

		pstn[0] = x
		pstn[1] = y 
		pstn[2] = z 
		
		for d0 in range(3):
			for d1 in range(3):
				Ktsr[d0,d1] = A[d0,d1]/chi + 1./3.*g[d0,d1]*K/chi
			
# ----------------------------------------------------------
# -------------- MADM measure ------------------------------
# ----------------------------------------------------------

		for d0 in range(3):
			for d1 in range(3):
				for d2 in range(3):
					for d3 in range(3):
 						MADMInt += s[d0]/(16.*np.pi)*gu[d1,d3]*gu[d0,d2]*(1./np.sqrt(chi)*(dgdx[d2,d3,d1] - dgdx[d1,d3,d2]) - 1./chi**(3./2.)*(g[d2,d3]*dchidx[d1] - g[d1,d3]*dchidx[d2]))
		
				
		MADM += w[k]*MADMInt*Rad*Rad*4.0*np.pi

		print("ADM mass = ", MADM)


# ----------------------------------------------------------
# -------------- JADM measure ------------------------------
# ----------------------------------------------------------
	
		for d0 in range(3):
			for d1 in range(3):
				for d2 in range(3):
					for d3 in range(3):
						JADMInt[d0] += - 1./(8*np.pi)*levi[d0,d1,d2]*s[d3]*pstn[d1]*(K*deltaij[d3,d2])
						for d4 in range(3):
							for d5 in range(3):
								JADMInt[d0] += + 1./(8*np.pi)*levi[d0,d1,d2]*s[d3]*pstn[d1]*(chi*chi*gu[d4,d3]*gu[d5,d2]*Ktsr[d4,d5]) 
			JADM[d0] += JADMInt[d0]*w[k]*Rad*Rad*4.0*np.pi

	print("MADM mass ", MADM, "at time",i.current_time)

	MADM_data.append(MADM)
	J1ADM_data.append(JADM[0])
	J2ADM_data.append(JADM[1])
	J3ADM_data.append(JADM[2])
	time_data.append(i.current_time)	
	CycleData.append(time.time()-start)
	np.savetxt('Cycletime.out',CycleData)
	np.savetxt('time.out',time_data)
	np.savetxt('MADM.out',MADM_data)
	np.savetxt('J1ADM.out',J1ADM_data)
	np.savetxt('J2ADM.out',J2ADM_data)
	np.savetxt('J3ADM.out',J3ADM_data)
