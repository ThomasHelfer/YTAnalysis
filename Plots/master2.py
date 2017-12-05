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
line = 2.5 # Line width
alp = 0.6 # alpha

#Loading dataset
ds = yt.load('../*.hdf5')

def checkneg(a):
        for i in np.arange(len(a)):
                if a[i]<0 :
                        return i

        return 0

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])

def _V(field, data):
	return 0.5*data["phi"]**2

def _Kinetic(field, data):
	return 0.5*data["Pi"]**2

def _Gradient(field, data):
	return 0.5*data["phi_gradient_magnitude"]**2

def _H2(field, data):
	return data["Ham"]*data["Ham"]

Gradient_data = []
Potential_data = []
Kinetic_data = []
Rho_data = []
Efolds_data = []
Gradient_data = []
L2H_data = []
time_data = []
CycleData = []

counter = 0

os.system("mkdir RhoPic")
os.system("mkdir ChiPic")
os.system("mkdir LapsePic")

for i in ds:

	counter += 1

        start = time.time()

	i.add_gradient_fields(('chombo', u'phi'))
	i.add_field("V",_V, units = "")
	i.add_field("Kinetic",_Kinetic, units = "")
	i.add_field("Gradient",_Gradient, units = "cm**(-2)")
	i.add_field("H2",_H2, units = "")

	for j in sorted(i.field_list):
  		print(j)

	ad = i.all_data()

	Rho_max = ad["rho"].max()

	time_data.append(i.current_time)

#==============================================
#       Normal plot - Rho
#==============================================


	Debug = time.time()

	plot_var = "rho"

        slc = yt.SlicePlot(i,'z',plot_var)
	slc.set_log(plot_var,False)
	slc.set_buff_size(1024)
#	slc.annotate_contour("chi")
        slc.set_cmap(field=plot_var, cmap='dusk')
	slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
	slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
	slc.set_colorbar_label("rho", r'$\rho \ \left[\frac{M_{pl}^2}{m}\right]$')
	slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
	#slc.set_window_size(10)
        slc.save()

	os.system("mv *.png RhoPic")

	print ("===============================================")
	print ("Time to compute RhoPlot ")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("===============================================")


#==============================================
#       Normal plot - Chi
#==============================================


	Debug = time.time()

	plot_var = "chi"

        slc = yt.SlicePlot(i,'z',plot_var)
	slc.set_log(plot_var,False)
	slc.set_buff_size(1024)
        slc.set_cmap(field=plot_var, cmap='dusk')
	slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
	slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
	slc.set_colorbar_label(plot_var, r'$\chi$')
	slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
	#slc.set_window_size(10)
        slc.save()

	os.system("mv *.png ChiPic")

	print ("===============================================")
	print ("Time to compute ChiPlot ")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("===============================================")


#==============================================
#       Normal plot - Lapse
#==============================================


	Debug = time.time()

	plot_var = "lapse"

        slc = yt.SlicePlot(i,'z',plot_var)
	slc.set_log(plot_var,False)
	slc.set_buff_size(1024)
        slc.set_cmap(field=plot_var, cmap='dusk')
	slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
	slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
	slc.set_colorbar_label(plot_var, r'$\alpha$')
	slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
	#slc.set_window_size(10)
        slc.save()

	os.system("mv *.png LapsePic")

	print ("===============================================")
	print ("Time to compute ChiPlot ")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("===============================================")
