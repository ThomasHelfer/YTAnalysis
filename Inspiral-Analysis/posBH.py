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
ds = yt.load('../outMatterSFp_*.hdf5')

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])
t0 = ds[0].current_time

time_data = []
CycleData = []
xData = []
yData = []
xData2 = []
yData2 = []
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

	c = i.r[center:,:,center]

	x = c["x"]
	y = c["y"]
	rho = c["chi"]	

	print (x)

	rho_min = float(rho.min())
	print (rho_min)
	x_max = np.where(rho == rho_min)[0][0]
	x_val = float(x[x_max])-center
	y_val = float(y[x_max])-center

	print("xval = ",x_val)
	print("yval = ",y_val)

	xData.append(x_val);
	yData.append(y_val);

	# -------------------------------
	#	line 
	# -------------------------------

	c = i.r[:center,:,center]

	x = c["x"]
	y = c["y"]
	rho = c["chi"]	

	print (x)

	rho_min = float(rho.min())
	print (rho_min)
	x_max = np.where(rho == rho_min)[0][0]
	x_val = float(x[x_max])-center
	y_val = float(y[x_max])-center

	print("xval = ",x_val)
	print("yval = ",y_val)

	xData2.append(x_val);
	yData2.append(y_val);
	
	
	time_data.append(i.current_time)
 	CycleData.append(time.time()-start)



N = len(xData)

xData = np.array(xData)
yData = np.array(yData)
xData2 = np.array(xData2)
yData2 = np.array(yData2)
time_data = np.array(time_data)

DT = time_data[1]-time_data[0]

xDeriv = (xData[:N-1]-xData[1:])/DT
yDeriv = (yData[:N-1]-yData[1:])/DT
xDeriv2 = (xData2[:N-1]-xData2[1:])/DT
yDeriv2 = (yData2[:N-1]-yData2[1:])/DT


rdata = np.sqrt(xData*xData+yData*yData)
rdata2 = np.sqrt(xData2*xData2+yData2*yData2)


Drdata = np.sqrt(xDeriv*xDeriv+yDeriv*yDeriv)
Drdata2 = np.sqrt(xDeriv2*xDeriv2+yDeriv2*yDeriv2)

#=====================================================
# Plots 
#=====================================================


plt.figure(figsize=(10,10))

thetaDat = np.linspace(0,2*3.15,100)
plt.plot(rdata[0]*np.sin(thetaDat),rdata[0]*np.cos(thetaDat),label = "Ref")

plt.scatter(yData,xData, alpha = alp)
plt.scatter(yData2,xData2, alpha = alp)
#plt.legend(loc = "lower right")
plt.xlabel(r'$x~[1/m]$')	
plt.ylabel(r'$y~[1/m]$')

plt.savefig("Trajec.pdf",bbox_inches = 'tight')
plt.close()	

plt.figure(figsize=(10,10))
plt.scatter(time_data,rdata, alpha = alp,label = "B1")
plt.scatter(time_data,rdata2, alpha = alp,label = "B2")
plt.legend(loc = "lower right")
plt.xlabel(r'$t~[1/m]$')	
plt.ylabel(r'$r~[1/m]$')

plt.savefig("TrajecR.pdf",bbox_inches = 'tight')
plt.close()	



plt.figure(figsize=(10,10))
plt.scatter(time_data[:N-1],Drdata, alpha = alp,label = "B1")
plt.legend(loc = "lower right")
plt.ylim([0,0.2])
plt.xlabel(r'$t~[1/m]$')	
plt.ylabel(r'$Dr~[1/m]$')

plt.savefig("TrajecDR.pdf",bbox_inches = 'tight')
plt.close()	



	
#	np.savetxt('Size.out',SizeData)
#	np.savetxt('time.out',time_data)
#	np.savetxt('Cycle.out',CycleData)
#	np.savetxt('maxRho.out',MasterData)
