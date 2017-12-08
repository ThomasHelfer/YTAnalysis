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

#yt.enable_parallelism()

# Plot parameters
line = 3 # Line width
alp = 0.6 # alpha
Quality = 'cubic' 

#Loading dataset
ds = yt.load('../GRChombo/GRChombo_1.0/S1x*.hdf5')

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])
t0 = ds[0].current_time

time_data = []
CycleData = []
MasterData = []
SizeData = []
counter = 0

for i in ds:

	counter += 1

        start = time.time()
	for j in sorted(i.field_list):
  		print(j)

	ad = i.all_data()

	# -------------------------------
	#	line 
	# -------------------------------

	c = i.ray((0,256,256),(256,256,256))

	x = c["x"]
	rho = c["rho"]	

	rho_max = float(rho.max())
	x_max = np.where(rho == rho_max)[0][0]
	rho_ren = rho-0.05*rho_max
	x_max_val = float(x[x_max])
	rhofun = interp1d(x, rho_ren, kind=Quality)
	x_1 = fsolve(rhofun,float(x[x_max])-1)[0]	
	x_2 = fsolve(rhofun,float(x[x_max])+1)[0]	

	size = x_2-x_1

	# -------------------------------
	#	slice  
	# -------------------------------

	sl = i.slice(2,256)
	rhosl = sl["rho"]
	xsl = sl["x"] 
	ysl = sl["y"] 
	rhosl_ren = rhosl - 0.05*rho_max
	index = np.where(rhosl_ren>0)
	
	xdata = xsl[index]
	ydata = ysl[index]

	index = np.where(xdata<256)

	xdata = xdata[index]
	ydata = ydata[index]
	
	MasterData.append(rho_max)
	time_data.append(i.current_time)
 	CycleData.append(time.time()-start)
 	SizeData.append(size)


	plt.figure(figsize=(20,20))
#	ax1 = plt.subplot(311)
#	plt.scatter(x,rho, alpha = alp,label = "t ="+str(i.current_time-t0))
#	plt.xlim([x_1-10,x_2+10])
#	plt.ylim([-0.0005,0.008])
#	plt.legend(loc = "upper right")
#	plt.xlabel(r'$t~[1/m]$')	
#	plt.ylabel(r'$\rho~[M_{pl}^2 m^2]$')

#	ax2 = plt.subplot(312)
	plt.scatter(xdata,ydata, alpha = alp,label = "t ="+str(i.current_time-t0))
#	plt.xlim([0,400])
	plt.xlim([x_max_val-10,x_max_val+10])
	plt.ylim([256-10,256+10])
	plt.legend(loc = "lower right")
	plt.xlabel(r'$x~[1/m]$')	
	plt.ylabel(r'$y~[1/m]$')

#	ax3 = plt.subplot(313)
#	plt.plot(time_data-t0,MasterData,linewidth = line,color = "darkviolet", alpha = alp,label = "rho_max")
#	plt.xlim([0,400])
#	plt.legend(loc = "lower right")
#	plt.xlabel(r'$t~[1/m]$')	
#	plt.ylabel(r'$\rho~[M_pl^2 m^2]$')
	plt.savefig(("rho%05d.png" % counter),bbox_inches = 'tight')
	plt.close()	
	
	np.savetxt('Size.out',SizeData)
	np.savetxt('time.out',time_data)
	np.savetxt('Cycle.out',CycleData)
	np.savetxt('maxRho.out',MasterData)
