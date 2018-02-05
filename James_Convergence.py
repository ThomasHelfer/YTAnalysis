# Convergence test
# James Cook

# Load the Modules
import sys
import yt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator
from yt import derived_field
import time
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Enable YT parallelism
yt.enable_parallelism()

# Plot Parameters
line = 2.5 # Line width
alp = 0.6 # alpha

########################
# Change Me
var = 'rho' # Variable to convergence test
Ord = 5
########################

####################################################
# Set-Up
####################################################

#mkdir the output directories
if not os.path.exists('data_convergence'):
    os.mkdir('data_convergence')

if not os.path.exists('data_plots'):
    os.mkdir('data_plots')

# Loading Datasets
# Low Res
low = yt.load("Low/GRChombo/Examples/Axion_Formation/*.3d.hdf5")

# Medium Res
medium = yt.load("Medium/GRChombo/Examples/Axion_Formation/*.3d.hdf5")

# High Res
high = yt.load("High/GRChombo/Examples/Axion_Formation/*.3d.hdf5")

# Set Up Paramters of Simulations
Simulation_Length = int(low[0].domain_right_edge[0])

Simulation_Low_Coarse = int(low[0].domain_dimensions[0])
Simulation_Medium_Coarse = 2. * Simulation_Low_Coarse
Simulation_High_Coarse = 2. * Simulation_Medium_Coarse

Simulation_Low_dx = Simulation_Length/float(Simulation_Low_Coarse)
Simulation_Medium_dx = Simulation_Length/float(Simulation_Medium_Coarse)
Simulation_High_dx = Simulation_Length/float(Simulation_High_Coarse)

low_centerXYZ =  low[0].domain_right_edge/2
low_center = float(low_centerXYZ[0])

medium_centerXYZ =  medium[0].domain_right_edge/2
medium_center = float(medium_centerXYZ[0])

high_centerXYZ =  high[0].domain_right_edge/2
high_center = float(high_centerXYZ[0])

# Arrays to keep data in
time_data = []
var_array_low = []
var_array_medium =[]
var_array_high =[]
convergence_array = []
convergence_lin_array =[]

# Dimension of interpolation array
dim = np.array([10,10,10])

# Central Point
central_point = medium_center + Simulation_Medium_dx/2.

# Some Basic Error Checking
# Check Grid Centers
if high_center != medium_center:
    print 'Grids have different centers!'
    print 'Aborting'
    sys.exit()
if high_center != low_center:
    print 'Grids have different centers!'
    print 'Aborting'
    sys.exit()
if medium_center != low_center:
    print 'Grids have different centers!'
    print 'Aborting'
    sys.exit()

# Check First Frames are Synchronised
if float(high[0].current_time) - float(medium[0].current_time) > 10e-15:
    print 'First frames are do not have the same time!'
    print 'Aborting'
    sys.exit()
if float(high[0].current_time) - float(low[0].current_time) > 10e-15:
    print 'First frames are do not have the same time!'
    print 'Aborting'
    sys.exit()
if float(medium[0].current_time) - float(low[0].current_time) > 10e-15:
    print 'First frames are do not have the same time!'
    print 'Aborting'
    sys.exit()

# Check amount of input frames
if len(high) != len(medium):
    print 'Different amount of frames!'
    print 'Aborting'
    sys.exit()
if len(high) != len(low):
    print 'Different amount of frames!'
    print 'Aborting'
    sys.exit()
if len(medium) != len(low):
    print 'Different amount of frames!'
    print 'Aborting'
    sys.exit()

####################################################
# Convergence Test
####################################################
for i in range(0,len(low)):
    # Check Frames are Synchronised
    if float(high[i].current_time) - float(medium[i].current_time) > 10e-5:
        print 'HM Frame ' + str(i) + ' are do not have the same time!'
        print high[i].current_time
        print medium[i].current_time
        print 'Aborting'
        sys.exit()
    if float(high[i].current_time) - float(low[i].current_time) > 10e-5:
        print 'HL Frame ' + str(i) + ' are do not have the same time!'
        print 'Aborting'
        sys.exit()
    if float(medium[i].current_time) - float(low[i].current_time) > 10e-5:
        print 'ML Frame ' + str(i) + ' are do not have the same time!'
        print 'Aborting'
        sys.exit()

    time_data.append(float(low[i].current_time))

    current_low_data = low[i]
    le = low_center - 5*Simulation_Low_dx
    low_covering_grid = current_low_data.covering_grid(level=0, left_edge=[le,le,le],dims=dim)

    xdata = np.array(low_covering_grid['x'][:,0,0])
    ydata = np.array(low_covering_grid['y'][0,:,0])
    zdata = np.array(low_covering_grid['z'][0,0,:])
    print xdata

    x0 = xdata[0]
    y0 = ydata[0]
    z0 = zdata[0]

    xcoord = (central_point-x0)/Simulation_Low_dx
    ycoord = (central_point-y0)/Simulation_Low_dx
    zcoord = (central_point-z0)/Simulation_Low_dx

    data = np.array(low_covering_grid[var])
    var_low = ndimage.map_coordinates(data, [[xcoord], [ycoord], [zcoord]], order=Ord)[0]
    var_array_low.append(var_low)

    interpolating_function = RegularGridInterpolator((xdata, ydata, zdata), data)
    var_low_lin = (interpolating_function([central_point,central_point,central_point])[0])

    low_covering_grid.clear_data()

    current_medium_data = medium[i]
    le = medium_center - 5*Simulation_Medium_dx
    medium_covering_grid = current_medium_data.covering_grid(level=0, left_edge=[le,le,le],dims=dim)

    xdata = np.array(medium_covering_grid['x'][:,0,0])
    ydata = np.array(medium_covering_grid['y'][0,:,0])
    zdata = np.array(medium_covering_grid['z'][0,0,:])
    print xdata

    x0 = xdata[0]
    y0 = ydata[0]
    z0 = zdata[0]

    xcoord = (central_point-x0)/Simulation_Medium_dx
    ycoord = (central_point-y0)/Simulation_Medium_dx
    zcoord = (central_point-z0)/Simulation_Medium_dx

    print xcoord

    data = np.array(medium_covering_grid[var])
    var_medium = ndimage.map_coordinates(data, [[xcoord], [ycoord], [zcoord]], order=Ord)[0]
    var_array_medium.append(var_medium)

    interpolating_function = RegularGridInterpolator((xdata, ydata, zdata), data)
    var_medium_lin = (interpolating_function([central_point,central_point,central_point])[0])

    medium_covering_grid.clear_data()

    current_high_data = high[i]
    le = high_center - 5*Simulation_High_dx
    high_covering_grid = current_high_data.covering_grid(level=0, left_edge=[le,le,le],dims=dim)

    xdata = np.array(high_covering_grid['x'][:,0,0])
    ydata = np.array(high_covering_grid['y'][0,:,0])
    zdata = np.array(high_covering_grid['z'][0,0,:])
    print xdata

    x0 = xdata[0]
    y0 = ydata[0]
    z0 = zdata[0]

    xcoord = (central_point-x0)/Simulation_High_dx
    ycoord = (central_point-y0)/Simulation_High_dx
    zcoord = (central_point-z0)/Simulation_High_dx

    print xcoord

    data = np.array(high_covering_grid[var])
    var_high = ndimage.map_coordinates(data, [[xcoord], [ycoord], [zcoord]], order=Ord)[0]
    var_array_high.append(var_high)

    interpolating_function = RegularGridInterpolator((xdata, ydata, zdata), data)
    var_high_lin = (interpolating_function([central_point,central_point,central_point])[0])

    high_covering_grid.clear_data()

    Convergence =  np.log2(abs(var_low-var_medium)/abs(var_medium-var_high))
    Convergence_lin =  np.log2(abs(var_low_lin-var_medium_lin)/abs(var_medium_lin-var_high_lin))
    convergence_array.append(Convergence)
    convergence_lin_array.append(Convergence_lin)

    # Write out
    np.savetxt('data_convergence/time.out',time_data)
    np.savetxt('data_convergence/var_high.out',var_array_high)
    np.savetxt('data_convergence/var_medium.out',var_array_medium)
    np.savetxt('data_convergence/var_low.out',var_array_low)
    np.savetxt('data_convergence/convergence_high.out',convergence_array)
    np.savetxt('data_convergence/convergence_lin.out',convergence_lin_array)

plt.plot(time_data,convergence_lin_array,linewidth = line, alpha = alp,color = 'indigo',label='Linear')
plt.ylim([0,4.3])
plt.ylabel('Order of Convergence')
plt.xlabel(r'$t\ [1/m_a]$')
plt.legend(loc = 'lower right')
#plt.legend()
plt.savefig("data_plots/Convergence_lin.png",bbox_inches = 'tight')
