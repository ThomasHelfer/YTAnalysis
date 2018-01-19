import yt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
from matplotlib import cm
import time
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

yt.enable_parallelism()

# Plot parameters
line = 2.5 # Line width
alp = 0.6 # alpha

#Loading dataset
ds = yt.load('../GRChombo/GRChombo_1.0/M1x*.hdf5')

def checkneg(a):
        for i in np.arange(len(a)):
                if a[i]<0 :
                        return i

        return 0

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])
max_extent = 20
counter = 0

os.system("mkdir SurfRend1")
os.system("mkdir VRStuff")

for i in ds:
	counter += 1

        start = time.time()
	ad = i.all_data()
	Rho_max = ad["rho"].max()
	C_max = ad["chi"].max()
	C_min = ad["chi"].min()


#==============================================
#	Surface_1 rendering Rho_max/2
#==============================================
	
	print("plotting.",i)

	Debug = time.time()

	#Isosurface value
	ISO = Rho_max/2.

        surface = i.surface(ad, "rho", ISO)
	print("plotting ...")
	colors = yt.apply_colormap(surface["chi"], cmap_name="hot")
	fig = plt.figure(figsize=(12,10))
	ax = fig.gca(projection='3d')
	print("plotting ...")
	p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)
	p3dc.set_facecolors(colors[0,:,:]/255.)
	ax.add_collection(p3dc)


	print("plotting ...")
#	max_extent = (surface.vertices.max(axis=1) - surface.vertices.min(axis=1)).max()
	centers = (surface.vertices.max(axis=1) + surface.vertices.min(axis=1)) / 2
	bounds = np.zeros([3,2])
	bounds[:,0] = centers[:] - max_extent/2
	bounds[:,1] = centers[:] + max_extent/2
	ax.auto_scale_xyz(bounds[0,:], bounds[1,:], bounds[2,:])


	m = cm.ScalarMappable(cmap="hot")
	m.set_array([C_min,C_max])
	m.set_clim(vmin=C_min,vmax=C_max)
	cbar = fig.colorbar(m, shrink=.5)
	cbar.set_label(r'$\chi$')

	ax.set_xlabel('x [1/m]')
	ax.set_ylabel('y [1/m]')
	ax.set_zlabel('z [1/m]')
	ax.set_title(r'Isocontour of max($\rho$)/2 = '+str(float(ISO)/2.)+ r' $\frac{M_{pl}^2}{m}$' )

	
	surface.export_obj("%s" % i , color_field='chi')
	os.system("mv *.mtl VRStuff")
	os.system("mv *.obj VRStuff")

	plt.savefig("%s_Surface.png" % i, dpi = 300)
	os.system("mv *.png SurfRend1")
	


	#surface.export_obj("%s_VR" % i )
	#os.system("mv *.obj *.mtl OBJRend1")

	print ("===============================================")
	print ("Time to compute SurfaceRender 1 ")
	print ("time = ",  time.time()-Debug ,"s ")
	print ("===============================================")

