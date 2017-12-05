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

# Enable YT parallelism
yt.enable_parallelism()

# Plot parameters
line = 3 # Line width
alp = 0.6 # alpha
Quality = 'cubic'

# Loading dataset
ds = yt.load('../*.hdf5')

#mkdir the output directories
if not os.path.exists('star_data'):
    os.mkdir('star_data')

if not os.path.exists('star_plot'):
    os.mkdir('star_plot')

# Find the center of the input data
centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])

# Define some useful functions
def _RhoJam(field, data):
        return data["rho"]/(data["chi"] * np.sqrt(data["chi"]))

# Define some arrays for the output from the calculation
time_data = []
star_size_data =[]
star_size_adjusted_data = []
rho_data =[]
rhojam_data = []

#Counter for output file names
counter=0

# Start the loop over frames
for i in ds:

    counter+=1
    
    # Add the RhoJam Field to the data
    i.add_field("RhoJam",_RhoJam, units = "")

    # Get the frame time
    time_data.append(i.current_time)

    # Line out from 0 to max (x), at the center of (y and z)
    c = i.ray((0,center,center),(center*2.0,center,center))

    # Name x coordinate and rho from lineout
    x = c["x"]
    rho = c["rho"]

    # Find maximum rho value within lineout
    rho_max = float(rho.max())

    # Find the coordinate of where the maximum value is
    x_max = np.where(rho == rho_max)[0][0]
    x_max_val = x[x_max]

    # Create a function that goes to negative when we are below 95% of rho max
    rhofun = interp1d(x, rho-0.05*rho_max, kind=Quality)

    # Solve the function for where it goes to zero using a best guess past and
    # before the the max value of x
    x_1 = fsolve(rhofun,float(x[x_max])-1)[0]
    x_2 = fsolve(rhofun,float(x[x_max])+1)[0]

    # Size of the 95% rho in center
    size = x_2-x_1
    star_size_data.append(size)

    # Create a sphere (in the full data) that has the same radius of the star
    sp = i.sphere([center, center, center], size/2.0)

    # Calculate mean Efolds, RhoTotal, and RhoJamTotal (integrals have been 'hacked')
    Efolds = float(sp.mean("Efolds", weight="cell_volume"))
    RhoTotal = float(sp.mean("rho", weight="cell_volume"))*(4./3.)*np.pi*np.power((size/2.),3)
    rho_data.append(RhoTotal)
    RhoJamTotal = float(sp.mean("RhoJam", weight="cell_volume"))*(4./3.)*np.pi*np.power((size/2.),3)
    rhojam_data.append(RhoJamTotal)

    # Calculate star size adjust for expansion
    star_size_adjusted_data.append(size*np.exp(Efolds))

    # Write out
    np.savetxt('star_data/time.out',time_data)
    np.savetxt('star_data/star_size.out',star_size_data)
    np.savetxt('star_data/star_size_adjusted.out',star_size_adjusted_data)
    np.savetxt('star_data/rho.out',rho_data)
    np.savetxt('star_data/rhojam.out',rhojam_data)

    plt.figure(figsize=(20,20))
    plt.scatter(x,rho, alpha = alp,label = "t ="+str(i.current_time))
    plt.xlim([x_1-5,x_2+5])
    plt.ylim([0,rho_max])
    plt.legend(loc = "upper right")
    plt.xlabel(r'$x~[1/m]$')
    plt.ylabel(r'$\rho~[M_{pl}^2 m^2]$')
    plt.savefig(("star_plot/rhojam%05d.png" % counter),bbox_inches = 'tight')
    plt.close()
