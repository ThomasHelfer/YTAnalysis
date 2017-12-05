import yt
#import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

#yt.enable_parallelism()


# ===============================================
#	 	 BH Horizon finder
# 	- assumes in spherical symmetry - 
#       Output : Black Hole mass evolution
# ===============================================

# Data file
loc = '../GRChombo/GRChombo_1.0/M1x00*'


def checkneg(a):
	for i in np.arange(len(a)):
		if a[i]<0 :
			return i
		
	return 0


ds = yt.load(loc)

# Set centerpoint of x - coord
centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])
# Inital guess for radius 
BHradguess = 1.5
# Quality = 'linear' , 'cubic' 
Quality = 'cubic' 

 
# Loading dataset

AHmark = True 
AHradii = []
AHradiitime = []
MassBHlist = []
MassBHlisty = []
Errorlist = []
# For error estimate 


for i in ds : 
	
	oray = i.ortho_ray(0, (float(centerXYZ[1]),float(centerXYZ[2])))
	c_Omega = oray["Omega"]
	c_x = oray["x"]

	i.print_stats()
	AHguess = checkneg(c_Omega)

	BHradguess = abs(float(c_x[AHguess])-center)

	if AHguess == 0 : 
		AHChecker = False
	else : 
		AHChecker = True

# Set time where BH Horizon appears 	
	if (AHmark and (AHChecker)):
		AHapperance = i.current_time
		AHmark = False

# Horizon finder 
	if ((AHChecker)): 	
		# Interpolating all needed functions 	
		OmegaInterpol = interp1d(c_x, c_Omega, kind=Quality)
		gamma22Interpol = interp1d(c_x, oray["gamma22"], kind=Quality)
		gamma33Interpol = interp1d(c_x, oray["gamma33"], kind=Quality)
		chiInterpol = interp1d(c_x, oray["chi"], kind=Quality)

		# Esimating resolution near horizon 
		Smallestres = float(c_x[AHguess]-c_x[AHguess-1])

		# Finding Apparent Horizon
		AHrad = fsolve(OmegaInterpol,center+BHradguess)	
		AHradii.append(AHrad-center)

		print ("BHradguess", BHradguess)
		print ("AHrad",AHrad -center )
		print ("BHradguess2", abs(float(c_x[AHguess-1])-center))

		# MassBH
		gamma22 = gamma22Interpol(AHrad)
		gamma33 = gamma33Interpol(AHrad)
		chi	= chiInterpol(AHrad)
		r0 = AHrad-center
		MassBH = 1./(2.*chi)*(r0)*(gamma22*gamma33)**(1./4.)
		MassBHlist.append(MassBH)

	
		Error = (1./(2.*chi)*(0.5*Smallestres)*(gamma22*gamma33)**(1./4.))

		# ==================================================
		# Calulating Error due to lack of spherical symmetry
		# ==================================================
		oray = i.ortho_ray(1, (float(centerXYZ[1]),float(centerXYZ[2])))
		c_Omega = oray["Omega"]
		c_y = oray["y"]
		OmegaInterpol = interp1d(c_y, c_Omega, kind=Quality)
		gamma22Interpol = interp1d(c_y, oray["gamma22"], kind=Quality)
		gamma33Interpol = interp1d(c_y, oray["gamma33"], kind=Quality)
		chiInterpol = interp1d(c_y, oray["chi"], kind=Quality)
                AHrad = fsolve(OmegaInterpol,center+BHradguess)
		gamma22 = gamma22Interpol(AHrad)
		gamma33 = gamma33Interpol(AHrad)
		chi	= chiInterpol(AHrad)
		r0 = AHrad-center
		MassBHy = 1./(2.*chi)*(r0)*(gamma22*gamma33)**(1./4.)
		MassBHlisty.append(MassBHy)

		print ("y-value mass ", MassBHy)

		Error += abs(MassBHy-MassBH)

		# ===================================================

		Errorlist.append(Error)
		# Current time 
		AHradiitime.append(i.current_time)


		np.savetxt('MassBHlist.out',MassBHlist)
		np.savetxt('MassBHlisty.out',MassBHlisty)
		np.savetxt('AHradiitime.out',AHradiitime)
		np.savetxt('Errorlist.out',Errorlist)


# Rescale to reduced Plank units (For Axion Star Paper)
#AHradiitime = np.array(AHradiitime)
#MassBHlist = np.array(MassBHlist)[:,0]
#Errorlist = np.array(Errorlist)[:,0]

#plk = np.sqrt(8*np.pi)
#AHradiitime = (1.0/plk)*AHradiitime
#MassBHlist =  plk*MassBHlist
#Errorlist = plk*Errorlist

# Horizon time output

print('---------------')
print('---------------')
print('Apparent horizon formed at time t = ', float(AHapperance),' 1/m ' )
print('---------------')
print('---------------')
## Plotting 
#plt.plot(AHradiitime,MassBHlist,linewidth = line, alpha = alp,color = 'blueviolet',label='Black hole mass')
#plt.fill_between(AHradiitime,(MassBHlist-Errorlist),(MassBHlist+Errorlist),linewidth = line, alpha = alp,color = 'skyblue')
#plt.legend(loc = 'lower right')
#plt.ylabel(r'$M_{BH}\ [M_{pl}^2/m_a]$')
#plt.xlabel(r'$t\ [1/m_a]$')
#plt.xlim([min(AHradiitime),max(AHradiitime)])
#plt.savefig("BlackHoleMass.png",bbox_inches = 'tight')
#plt.show()
