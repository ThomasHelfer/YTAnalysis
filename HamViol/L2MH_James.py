# L2MH_James.py
# L2H and L2M calculater based of HamViol.py
# Different tab aligning to HamViol.py

# Load the modules
import yt
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Enable Parallelism
yt.enable_parallelism()

# mkdir the plot directory
if not os.path.exists('plots'):
    os.mkdir('plots')

# Plot parameters
line = 2.5 # Line width
alp = 0.6 # alpha

# CHANGE ME
# Loading dataset (Load from one folder up)
ds = yt.load('../outMatterSF_*.3d.hdf5')

# Get center of 0th frame of simulation
centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])

# Define H2
def _H2(field, data):
        return data["Ham"]*data["Ham"]

# Define M2
def _M2(field, data):
        return data["Mom1"]*data["Mom1"] + data["Mom2"]*data["Mom2"] + data["Mom3"]*data["Mom3"]

# Generate our arrays to store the results in
time_data = []
L2H_data = []
L2M_data = []
CycleData = []

for i in ds:
    start = time.time()

    i.add_field("H2",_H2, units = "")
    i.add_field("M2",_M2, units = "")

    ad = i.r[8:120,8:120,8:120]

    time_data.append(i.current_time)

    #########
    # L2H
    #########

    meanH2 = ad.mean("H2", weight="cell_volume")
    L2H = np.sqrt(meanH2)
    L2H_data.append(L2H)

    #########
    # L2M
    #########

    meanM2 = ad.mean("M2", weight="cell_volume")
    L2M = np.sqrt(meanM2)
    L2M_data.append(L2M)

    #########
    # Write Out
    #########

    CycleData.append(time.time()-start)

    np.savetxt('time.out',time_data)
    np.savetxt('cycle.out',CycleData)
    np.savetxt('L2M.out',L2M_data)
    np.savetxt('L2H.out',L2H_data)

# Plots
# Update rcParams so that we get all axis labelling
rcParams.update({'figure.autolayout': True})
rcParams['axes.formatter.limits'] = [-5,5]

# Graph Out
plt.figure(1)
plt.plot(time_data,L2H_data)
plt.title('$\\mathcal{H}$ vs time')
plt.ylabel('$\\mathcal{H}$')
plt.xlabel('Relax Time')
plt.savefig('plots/H.png')

# Graph Out
plt.figure(2)
plt.plot(time_data,L2M_data)
plt.title('$\\mathcal{M}$ vs time')
plt.ylabel('$\\mathcal{M}$')
plt.xlabel('Relax Time')
plt.savefig('plots/M.png')
