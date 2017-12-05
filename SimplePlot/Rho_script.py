import yt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
import time 
import os 
yt.enable_parallelism()

# ===============================================
#	    Pictures
# ===============================================

# Plot parameters 
line = 2.5 # Line width
alp = 0.6 # alpha


# EXTRA CODE FOR REDUCED PLANK UNITS
# -----------------------------------------------
# -----------------------------------------------
#units_override = {"length_unit":(1.0/plk,"l_pl"),
#                  "time_unit":(1.0/plk,"t_pl"),
#                  "mass_unit":(plk,"m_pl")}
#ds = yt.load(loc, units_override=units_override)
# -----------------------------------------------
# -----------------------------------------------

# Data file
loc = '../GRChombo/GRChombo_1.0/S1x*'


#Loading dataset

ds = yt.load(loc)

SimStarttime = ds[0].current_time 
cycledata = []
counter = 0 

for i in ds:

	start = time.time()
        slc = yt.SlicePlot(i,'z','rho',width = (90.0, 'cm'))
#	slc.annotate_grids()
	slc.set_log('rho',False)
	slc.set_buff_size(1024)
#	slc.set_zlim('rho', 0.0007, 0.0000)
        slc.set_cmap(field="rho", cmap='dusk')
	slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
	slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
	slc.set_colorbar_label("rho", r'$\rho \ \left[\frac{M_{pl}^2}{m}\right]$')
	slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time-SimStarttime))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
	slc.set_window_size(10)
        slc.save()

	cycletime = time.time()	- start 
	cycledata.append(cycletime)
	np.savetxt('cycledata.out',cycledata)
	
	counter += 1 	

	
