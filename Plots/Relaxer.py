# James Cook
# Relaxer.py
# Nov 2017
# A script to calculate optimum relax time for Axion Star Formation

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

program_time = time.time()

# Plot parameters
line = 2.5 # Line width
alp = 0.6 # alpha

# Useful Constants
eightpi = 8*np.pi

# Set up first output File
print 'Output Summary'
print ' '

# Load the Files
start_time = time.time()
ds = yt.load('../*.hdf5')
end_time = time.time() - start_time

print 'File Load: ' ,end_time, ' Seconds'

time_data = []
Viol_array = []
Rho_tot = []
Rho_data = []
K_data = []

# Begin File Analysis
for i in ds:

    ad = i.all_data()

    # Get Current Time
    time_data.append(i.current_time)

    # Get Max Rho
    Rho_max = ad["rho"].max()

    # Get Max Ham
    Ham_max = ad["Ham"].max()
    Ham_min = ad["Ham"].min()

    if np.sqrt(Ham_min*Ham_min)>Ham_max:
        Ham_max = Ham_min

    # Calculate Relative Violation
    Viol = Ham_max/(16*np.pi*Rho_max)*100
    Viol_array.append(Viol)

    # Total Rho
    Rho_tot.append(ad.quantities.total_quantity(["rho"]))

    # Mean Total Energy
    Rho = float(ad.mean("rho", weight="cell_volume"))
    Rho_data.append(Rho_max)

    Inflation_K = -(3.0*eightpi*Rho)**(0.5)
    K_data.append(Inflation_K)

#famper = interp1d(Viol_array, time_data, 'cubic', fill_value='extrapolate')
#fh = interp1d(time_data,Viol_array, 'cubic', fill_value='extrapolate')

# Write Out
np.savetxt('time.out',time_data)
np.savetxt('viol.out',Viol_array)
np.savetxt('inflation_K.out',K_data)
np.savetxt('rhomax.out',Rho_data)

#print 'Optimum Relax Time', famper(0.001)

#Smooth_time= np.linspace(0,time_data[-1],100)

# Graph Out
plt.figure(1)
plt.plot(time_data,Viol_array)
#plt.plot(Smooth_time,fh(Smooth_time))
plt.title('$\\mathcal{H}$ vs time')
plt.ylabel('$\\mathcal{H}$')
plt.xlabel('Relax Time')
plt.savefig('Viol.png')

print 'Program Time: ' ,time.time() - program_time, ' Seconds'
