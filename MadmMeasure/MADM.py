import yt
from yt import derived_field
import numpy as np
from numpy.linalg import inv
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import time
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
ds = yt.load('../../*.hdf5')

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])

RadiusADM = 100

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
BHmass = []
BHmassError = []

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

	c = i.point([1,2,3])
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


	print ("Calculating Theta ... ")
	MADM = np.zeros(N)

# ----------------------------------------------------------
	for i in range(3):
		for j in range(3):
			MADM += (dgdx[i,j,i]/chi-dchidx[i]/chi**2.*g[i,j])*s[j]


	for i in range(3):
		for j in range(3):
			MADM += (-chi**(-2.)*dchidx[i]*g[j,j]+dgdx[j,j,i]*chi**(-1.))*s[i]


	CycleData.append(time.time()-start)
	np.savetxt('Cycletime.out',CycleData)
