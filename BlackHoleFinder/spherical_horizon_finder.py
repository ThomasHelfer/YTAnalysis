import yt
from yt import derived_field
import numpy as np
from numpy.linalg import inv
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import time
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
yt.enable_parallelism()

# ==============================================================================
#		Spherical Horizon finder 
#		
# Assumes (Obviously) spherical symmetry !! Gives bad results if not !
# Input = Black hole center
# ==============================================================================

#Loading dataset
ds = yt.load('../out*.hdf5')

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])


# =================================
#  INPUT #  
# =================================
BHcenter = [center,center,center]
# ================================
# ================================

time_data = []
CycleData = []
AHradius = []
AHradiusError = []


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

	c = i.ray((BHcenter[0],BHcenter[1],BHcenter[2]),(BHcenter[0]*2,BHcenter[1],BHcenter[2]))
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
	gu = [[],[],[]]

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

	Omega = c["omega"][srt]

	Kscalar = c["K"][srt]
	Kscalar = np.array(Kscalar)	
	
	A11 = np.array(c["A11"][srt]);
	A12 = np.array(c["A12"][srt]);
	A13 = np.array(c["A13"][srt]);
	A22 = np.array(c["A22"][srt]);
	A23 = np.array(c["A23"][srt]);
	A33 = np.array(c["A33"][srt]);

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
	gu = np.zeros([3,3,N])
	dsdx = np.zeros([3,3,N])
	s = np.zeros([3,N])

	for IND in range(N):
		# Definitions 
		# ---------------------	
		s[0,IND] = x[IND]/r[IND]
		s[1,IND] = y[IND]/r[IND]
		s[2,IND] = z[IND]/r[IND] 
		dsdx[0,0,IND] = y[IND]*y[IND]/(pow(r[IND],3)) + z[IND]*z[IND]/(pow(r[IND],3)) 
		dsdx[1,1,IND] = x[IND]*x[IND]/(pow(r[IND],3)) + z[IND]*z[IND]/(pow(r[IND],3))
		dsdx[2,2,IND] = x[IND]*x[IND]/(pow(r[IND],3)) + y[IND]*y[IND]/(pow(r[IND],3))
		dsdx[0,1,IND] = -s[0,IND]*s[1,IND]/r[IND]
		dsdx[1,0,IND] = dsdx[0,1,IND]
		dsdx[0,2,IND] = -s[0,IND]*s[2,IND]/r[IND]
		dsdx[2,0,IND] = dsdx[0,2,IND]
		dsdx[1,2,IND] = -s[1,IND]*s[2,IND]/r[IND]
		dsdx[2,1,IND] = dsdx[1,2,IND]
 
		# Inversion  
		# ---------------------	
		temp = g[:,:,IND]
		tempINV = inv(temp)
		gu[:,:,IND] = tempINV

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


	print ("Calculating Theta ... ")
	E = np.zeros(N)
	B = np.zeros(N)
	C = np.zeros(N)
	D = np.zeros(N)
	Theta = np.zeros(N)

	for d in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					E += - chi*chi*gu[d,k]*s[k]*gu[j,l]*s[l]*dsdx[j,d] - 0.5*chi*gu[d,j]*s[j]*(gu[k,l]*dchidx[d])*s[k]*s[l];
					for n in range(3):
						for m in range(3):
							E += 0.5*chi*gu[d,j]*s[j]*chi*gu[k,m]*gu[l,n]*dgdx[m,n,d]*s[k]*s[l];

	for d in range(3):
		for j in range(3):
			B += chi*gu[d,j]*dsdx[j,d] - 3.0/2.0*dchidx[d]*gu[d,j]*s[j] + gu[d,j]*dchidx[d]*s[j]; 
			for k in range(3):
				for l in range(3):
					B += - chi*gu[d,k]*gu[j,l]*dgdx[k,l,d]*s[j]; 



	for d in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					C += chi*gu[d,k]*gu[j,l]*(A[k,l]+1.0/(3.0)*Kscalar*g[k,l])*s[d]*s[j];

	for d in range(3):
		for j in range(3):
			D += chi*gu[d,j]*s[d]*s[j];

	Theta = 1.0/np.sqrt(D)*(E/D + B) + C/D - Kscalar;
		
		
	

# =====================================================
#	Finding Horizon
# =====================================================
		
	AHguess = 0 
	for d in range(N):
		if Theta[N-1-d]<0:
			AHguess = N-1-d
			break

	if AHguess == 0 : 
		AHChecker = False
	else : 
		AHChecker = True
		print("Black hole horizon has formed !!")
		time_data.append(i.current_time)
		

	if (AHChecker):
		ThetaInterpol = interp1d(r, Omega, kind=Quality)
		# Estimating Error by using local resolution 
		AH_radius_error = abs(r[AHguess]-r[AHguess-1])
		AHradiusError.append(AH_radius_error)
		BHradguess = r[AHguess]
		AHrad = fsolve(ThetaInterpol,BHradguess+AH_radius_error)	
		print (BHradguess)
		AHradius = []


		
		plt.figure(figsize=(20,10)) 	
		plt.scatter(r,Theta,label = "Calculation",linewidth = 3,color = "darkviolet",linestyle = "--")
		plt.scatter(r,Omega,label = "Ref",linewidth = 3,color = "black",linestyle = "--")
		plt.plot(np.ones(N)*AHrad,np.linspace(min(Theta)*2,max(Theta)*2,N),label = "Black hole horizon",linewidth = 1.4,color = "black")
		plt.legend(loc = "best")
		plt.ylim([min(Theta)*2,max(Theta)*2])
		plt.ylabel(r"Quality $[\%]$")
		plt.xlabel(r"r$[M]$")
		plt.xlim([0,50])
		plt.savefig("ADR.png")
		plt.show()


	CycleData.append(time.time()-start)
	np.savetxt('Cycletime.out',CycleData)
	np.savetxt('AH_radius_error.out',AHradiusError)
	np.savetxt('black_hole_formation_time.out',time_data)
