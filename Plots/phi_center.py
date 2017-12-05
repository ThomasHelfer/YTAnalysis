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

yt.enable_parallelism()

# Plot parameters
line = 3 # Line width
alp = 0.6 # alpha
Quality = 'cubic'

#Loading dataset
ds = yt.load('../*.hdf5')

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])
t0 = ds[0].current_time

time_data = []
CycleData = []
MasterData = []


for i in ds:


        start = time.time()
	for j in sorted(i.field_list):
  		print(j)

	ptn = i.point(centerXYZ)
	phi = float(ptn["phi"])

	MasterData.append(phi)
	time_data.append(i.current_time)

	np.savetxt('time1.out',time_data)
	np.savetxt('phi.out',MasterData)
