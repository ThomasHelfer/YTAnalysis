import yt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
import time

#yt.enable_parallelism()

# Plot parameters 
line = 2.5 # Line width
alp = 0.6 # alpha

#Loading dataset

ds = yt.load('../GRChombo/GRChombo_1.0/M1x*.3d.hdf5')

def checkneg(a):
        for i in np.arange(len(a)):
                if a[i]<0 :
                        return i

        return 0

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])

@derived_field(name = "H2", units = "")
def _H2(field, data):
	return data["H"]*data["H"]

Madm_data = []
time_data = []
L2H_data = []
CycleData = []

for i in ds:

        start = time.time()

	ad = i.r[8:504,8:504,8:504]
	hot_ad = ad.cut_region(["obj['Omega'] > 0.0"])

        time_data.append(i.current_time)
#==============================================
#	Madm 
#==============================================
	Debug = time.time()
	vol_Madm = ad.sum("cell_volume")
	print ("voltime = ", Debug - time.time())
	Debug = time.time()
 	Madm = ad.sum("Madm")
	print ("voltime = ", Debug - time.time())
	Madm_data.append(Madm)

#==============================================
#	L2H
#==============================================

        meanH2 = hot_ad.mean("H2", weight="cell_volume")
	L2H = np.sqrt(meanH2)	
	L2H_data.append(L2H)


#==============================================
#	write out  
#==============================================
 	CycleData.append(time.time()-start)	

	np.savetxt('time.out',time_data)
	np.savetxt('Madm.out',Madm_data)
	np.savetxt('L2H.out',L2H_data)


