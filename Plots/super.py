# James Cook
# super.py
# Nov 2017
# A script to fully analyse Axion Star Formation Files

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
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

yt.enable_parallelism()

program_time = time.time()

# mkdir the plot directory
if not os.path.exists('plots'):
    os.mkdir('plots')

#mkdir the data directory
if not os.path.exists('data'):
    os.mkdir('data')

#mkdir the data directory
if not os.path.exists('RhoPic'):
    os.mkdir('RhoPic')

#mkdir the data directory
if not os.path.exists('ChiPic'):
    os.mkdir('ChiPic')

#mkdir the data directory
if not os.path.exists('LapsePic'):
    os.mkdir('LapsePic')

#mkdir the data directory
if not os.path.exists('HamPic'):
    os.mkdir('HamPic')

# Plot parameters
line = 2.5 # Line width
alp = 0.6 # alpha
Quality = 'cubic'

#Specify Relax Time
relax_time = 600 # CHANGE ME
counter = 0 # Counter to index data so that relax time can be removed

# Set up first output File
print 'Output Summary'
print ' '

# Load the Files
ds = yt.load('../*.hdf5')

centerXYZ =  ds[0].domain_right_edge/2
center = float(centerXYZ[0])
t0 = ds[0].current_time

def _V(field, data):
        return 0.5*data["phi"]**2

def _Kinetic(field, data):
        return 0.5*data["Pi"]**2

def _Gradient(field, data):
        return 0.5*data["phi_gradient_magnitude"]**2

def _H2(field, data):
        return data["Ham"]*data["Ham"]

# Create the Arrays for the Data
time_data = []
Viol_array = []
phi_center_data = []
Gradient_data = []
Potential_data = []
Kinetic_data = []
Rho_data = []
Rho_tot = []
Efolds_data = []
Gradient_data = []
L2H_data = []
CycleData = []

# Begin File Analysis
for i in ds:

    i.add_gradient_fields(('chombo', u'phi'))
    i.add_field("V",_V, units = "")
    i.add_field("Kinetic",_Kinetic, units = "")
    i.add_field("Gradient",_Gradient, units = "cm**(-2)")
    i.add_field("H2",_H2, units = "")

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

    # Calculate Phi at the Center
    ptn = i.point(centerXYZ)
    phi = float(ptn["phi"])
    phi_center_data.append(phi)

    # Mean Potential Energy Density
    Pot = ad.mean("V", weight="cell_volume")
    Potential_data.append(Pot)

    # Mean potential energydensity
    Kin = float(ad.mean("Kinetic", weight="cell_volume"))
    Kinetic_data.append(Kin)

    # Mean gradient energydensity
    Grad = float(ad.mean("Gradient", weight="cell_volume"))
    Gradient_data.append(Grad)

    # Mean Total Energy
    Rho = float(ad.mean("rho", weight="cell_volume"))
    Rho_data.append(Rho)

    # Mean Total L2H
    meanH2 = ad.mean("H2", weight="cell_volume")
    L2H = np.sqrt(meanH2)
    L2H_data.append(L2H)

    # Mean Total Efolds
    Efolds = float(ad.mean("Efolds", weight="cell_volume"))
    Efolds_data.append(Efolds)

    # Write Out
    np.savetxt('data/time.out',time_data)
    np.savetxt('data/HRelRelax.out',Viol_array)
    np.savetxt('data/phic.out',phi_center_data)
    np.savetxt('data/Potential.out',Potential_data)
    np.savetxt('data/Kinetic.out',Kinetic_data)
    np.savetxt('data/Gradient.out',Gradient_data)
    np.savetxt('data/Efolds.out',Efolds_data)
    np.savetxt('data/rho.out',Rho_data)
    np.savetxt('data/L2H.out',L2H_data)
    np.savetxt('data/rhotot.out',Rho_tot)

    # Normal plot - Rho
    plot_var = "rho"
    slc = yt.SlicePlot(i,'z',plot_var)
    slc.set_log(plot_var,False)
    slc.set_buff_size(1024)
    slc.set_cmap(field=plot_var, cmap='dusk')
    slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
    slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
    slc.set_colorbar_label("rho", r'$\rho \ \left[\frac{M_{pl}^2}{m}\right]$')
    slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
    slc.save()
    os.system("mv *.png RhoPic")

    # Normal plot - Chi
    plot_var = "chi"
    slc = yt.SlicePlot(i,'z',plot_var)
    slc.set_log(plot_var,False)
    slc.set_buff_size(1024)
    slc.set_cmap(field=plot_var, cmap='dusk')
    slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
    slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
    slc.set_colorbar_label(plot_var, r'$\chi$')
    slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
    slc.save()
    os.system("mv *.png ChiPic")

    # Normal plot - Lapse
    plot_var = "lapse"
    slc = yt.SlicePlot(i,'z',plot_var)
    slc.set_log(plot_var,False)
    slc.set_buff_size(1024)
    slc.set_cmap(field=plot_var, cmap='dusk')
    slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
    slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
    slc.set_colorbar_label(plot_var, r'$\alpha$')
    slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
    slc.save()
    os.system("mv *.png LapsePic")

    # Normal plot - Ham
    plot_var = "Ham"
    slc = yt.SlicePlot(i,'z',plot_var)
    slc.set_log(plot_var,False)
    slc.set_buff_size(1024)
    slc.set_cmap(field=plot_var, cmap='dusk')
    slc.set_xlabel(r'x $\left[\frac{1}{m}\right]$')
    slc.set_ylabel(r'y $\left[\frac{1}{m}\right]$')
    slc.set_colorbar_label(plot_var, r'Ham')
    slc.annotate_text((0.13, 0.92),( 'time = '+str(float(i.current_time))+' 1/m' ) , coord_system='figure',text_args={'color':'white'})
    slc.save()
    os.system("mv *.png HamPic")

################################################################################
# Post yt Analysis
################################################################################

# Import the data files
time_data = np.genfromtxt('data/time.out')
e_folds_data = np.genfromtxt('data/Efolds.out')
rho_data = np.genfromtxt('data/rho.out')
l2h_data = np.genfromtxt('data/L2H.out')
gradient_data = np.genfromtxt('data/Gradient.out')
kinetic_data = np.genfromtxt('data/Kinetic.out')
potential_data = np.genfromtxt('data/Potential.out')

# Find index of relax_time
for i in time_data:
    if i <=relax_time:
        counter +=1

print 'End of Relax Time HRel = ', Viol_array[counter], '%'

# New array with only post relax time
post_relax_time = []

# Modify post relax time to start from zero
for i in time_data[counter:]:
    post_relax_time.append(i-relax_time)

# New array to calculate scale factor from Efolds data
scale_factor = []

for i in e_folds_data:
    scale_factor.append(np.exp(i))

# Comparison data for plots
rho_matter = []
rho_radiation = []
pot_matter =[]
pot_radiation = []
pot_scalar = []
pot_late_matter =[]
pot_late_radiation = []
pot_late_scalar = []
f_late_scalar=[]
k_late_matter =[]
k_late_radiation = []
k_late_scalar = []
k_late_f = []

# Scaling for converging at the end of the plots
m_scaling = potential_data[-1]/((1./scale_factor[-1])**3)
r_scaling = potential_data[-1]/((1./scale_factor[-1])**4)
f_scaling = potential_data[-1]/((1./scale_factor[-1])**5)
s_scaling = potential_data[-1]/((1./scale_factor[-1])**6)

mk_scaling = kinetic_data[-1]/((1./scale_factor[-1])**3)
rk_scaling = kinetic_data[-1]/((1./scale_factor[-1])**4)
fk_scaling = kinetic_data[-1]/((1./scale_factor[-1])**5)
sk_scaling = kinetic_data[-1]/((1./scale_factor[-1])**6)

for i in scale_factor[counter:]:
    rho_matter.append(rho_data[counter] * (1./i)**3)
    rho_radiation.append(rho_data[counter] * (1./i)**4)
    pot_matter.append(potential_data[counter] * (1./i)**3)
    pot_late_matter.append(m_scaling * (1./i)**3)
    pot_radiation.append(potential_data[counter] * (1./i)**4)
    pot_late_radiation.append(r_scaling* (1./i)**4)
    pot_scalar.append(potential_data[counter] * (1./i)**6)
    pot_late_scalar.append(s_scaling * (1./i)**6)
    f_late_scalar.append(f_scaling * (1./i)**5)
    k_late_matter.append(mk_scaling * (1./i)**3)
    k_late_radiation.append(rk_scaling* (1./i)**4)
    k_late_scalar.append(sk_scaling * (1./i)**6)
    k_late_f.append(fk_scaling * (1./i)**5)

# Plots
# Update rcParams so that we get all axis labelling
rcParams.update({'figure.autolayout': True})
rcParams['axes.formatter.limits'] = [-5,5]

# Quality of V K and Grad
plt.figure(1)
plt.title('Confidence in approximations')
plt.plot(post_relax_time,rho_data[counter:]-(kinetic_data[counter:]+gradient_data[counter:]+potential_data[counter:]))
plt.xlabel('$t\ [1/m]$')
plt.ylabel("$\\rho - V - K - \\nabla$")
plt.savefig('plots/confidence.png')

# Efolds
plt.figure(2)
plt.title('Efolds')
plt.plot(post_relax_time,e_folds_data[counter:])
plt.xlabel('$t\ [1/m]$')
plt.ylabel('Efolds')
plt.savefig('plots/efolds.png')

# L2H
plt.figure(3)
plt.title('Change in $L^2$')
plt.plot(post_relax_time,l2h_data[counter:])
plt.xlabel('$t\ [1/m]$')
plt.ylabel('$L^2 (H)$')
plt.savefig('plots/l2norm.png')

# Energy distribution
plt.figure(4)
plt.title('Change in energy distributions')
plt.plot(post_relax_time,gradient_data[counter:], label="$\\nabla$")
plt.plot(post_relax_time,kinetic_data[counter:] , label = "$K$")
plt.plot(post_relax_time,potential_data[counter:],label = "$V$")
plt.plot(post_relax_time,kinetic_data[counter:]+gradient_data[counter:]+potential_data[counter:],label = "$Total$")
plt.xlabel('$t\ [1/m]$')
plt.legend()
plt.savefig('plots/energy.png')

# Energy distribution (normalised)
plt.figure(5)
plt.title('Change in normalised energy distributions')
plt.plot(post_relax_time,gradient_data[counter:]/(kinetic_data[counter:]+gradient_data[counter:]+potential_data[counter:]), label="$\\nabla$")
plt.plot(post_relax_time,kinetic_data[counter:]/(kinetic_data[counter:]+gradient_data[counter:]+potential_data[counter:]) , label = "$K$")
plt.plot(post_relax_time,potential_data[counter:]/(kinetic_data[counter:]+gradient_data[counter:]+potential_data[counter:]),label = "$V$")
plt.plot(post_relax_time,(kinetic_data[counter:]+gradient_data[counter:]+potential_data[counter:])/(kinetic_data[counter:]+gradient_data[counter:]+potential_data[counter:]),label = "$Total$")
plt.xlabel('$t\ [1/m]$')
plt.legend()
plt.savefig('plots/norm_energy.png')

# L2H over relax time
plt.figure(6)
plt.plot(time_data[:counter],l2h_data[:counter])
plt.title('Change in $L^2$ vs relaxation time')
plt.xlabel('$t\ [1/m]$')
plt.ylabel('$L^2 (H)$')
plt.savefig('plots/relax_l2norm.png')

# Scale Factor
plt.figure(7)
plt.plot(post_relax_time,scale_factor[counter:])
plt.xlabel('$t\ [1/m]$')
plt.ylabel('$a(t)$')
plt.savefig('plots/scale_factor.png')

# Energy Dist vs Scale Factor
plt.figure(8)
plt.title('Change in energy distributions')
plt.plot(np.log(scale_factor[counter:]),gradient_data[counter:], label="$\\nabla$")
plt.plot(np.log(scale_factor[counter:]),kinetic_data[counter:] , label = "$K$")
plt.plot(np.log(scale_factor[counter:]),potential_data[counter:],label = "$V$")
plt.plot(np.log(scale_factor[counter:]),kinetic_data[counter:]+gradient_data[counter:]+potential_data[counter:],label = "$Total$")
plt.xlabel('$ln(a(t))$')
plt.legend()
plt.savefig('plots/energy_scale.png')

# rho_total behaviour
plt.figure(9)
plt.plot(scale_factor[counter:],rho_data[counter:], label = '$\\rho$ data')
plt.plot(scale_factor[counter:],rho_matter, label = '$\\rho$ matter', linestyle = 'dashed' )
plt.plot(scale_factor[counter:],rho_radiation, label = '$\\rho$ radiation', linestyle = 'dashed' )
plt.ylabel('$\\rho(a(t))$')
plt.xlabel('$a(t)$')
plt.legend()
plt.savefig('plots/type.png')

# Kinetic
plt.figure(10)
plt.plot(scale_factor[counter:],kinetic_data[counter:])
plt.ylabel('$K(a(t))$')
plt.xlabel('$a(t)$')
plt.savefig('plots/Kinetic.png')

# Potential
plt.figure(11)
plt.plot(scale_factor[counter:],potential_data[counter:], label = '$V(a(t))$ data')
plt.plot(scale_factor[counter:],pot_matter, label = '$n=-3$', linestyle = 'dashed' )
plt.plot(scale_factor[counter:],pot_radiation, label = '$n=-4$', linestyle = 'dashed' )
plt.plot(scale_factor[counter:],pot_scalar, label = '$n=-6$', linestyle = 'dashed' )
plt.ylabel('$V(a(t))$')
plt.xlabel('$a(t)$')
plt.legend()
plt.savefig('plots/Potential.png')

# Potential scaled at end
plt.figure(12)
plt.plot(scale_factor[counter:],potential_data[counter:], label = '$V(a(t))$ data')
plt.plot(scale_factor[counter:],pot_late_matter, label = '$n=-3$', linestyle = 'dashed' )
plt.plot(scale_factor[counter:],pot_late_radiation, label = '$n=-4$', linestyle = 'dashed' )
plt.plot(scale_factor[counter:],f_late_scalar, label = '$n=-5$', linestyle = 'dashed' )
plt.plot(scale_factor[counter:],pot_late_scalar, label = '$n=-6$', linestyle = 'dashed' )
plt.ylabel('$V(a(t))$')
plt.xlabel('$a(t)$')
plt.legend()
plt.savefig('plots/Potential_2.png')

# Potential scaled at end
plt.figure(13)
plt.plot(scale_factor[counter:],kinetic_data[counter:], label = '$V(a(t))$ data')
plt.plot(scale_factor[counter:],k_late_matter, label = '$n=-3$', linestyle = 'dashed' )
plt.plot(scale_factor[counter:],k_late_radiation, label = '$n=-4$', linestyle = 'dashed' )
plt.plot(scale_factor[counter:],k_late_f, label = '$n=-5$', linestyle = 'dashed' )
plt.plot(scale_factor[counter:],k_late_scalar, label = '$n=-6$', linestyle = 'dashed' )
plt.ylabel('$K(a(t))$')
plt.xlabel('$a(t)$')
plt.legend()
plt.savefig('plots/Kinetic_2.png')

# Potential scaled at end
plt.figure(14)
plt.plot(np.log(scale_factor[counter:]),np.log(potential_data[counter:]), label = '$V(a(t))$ data')
plt.plot(np.log(scale_factor[counter:]),np.log(pot_late_matter), label = '$n=-3$', linestyle = 'dashed' )
plt.plot(np.log(scale_factor[counter:]),np.log(pot_late_radiation), label = '$n=-4$', linestyle = 'dashed' )
plt.plot(np.log(scale_factor[counter:]),np.log(f_late_scalar), label = '$n=-5$', linestyle = 'dashed' )
plt.plot(np.log(scale_factor[counter:]),np.log(pot_late_scalar), label = '$n=-6$', linestyle = 'dashed' )
plt.ylabel('$V(a(t))$')
plt.xlabel('$a(t)$')
plt.legend()
plt.savefig('plots/Potential_2_log.png')

# Potential scaled at end
plt.figure(15)
plt.plot(np.log(scale_factor[counter:]),np.log(kinetic_data[counter:]), label = '$V(a(t))$ data')
plt.plot(np.log(scale_factor[counter:]),np.log(k_late_matter), label = '$n=-3$', linestyle = 'dashed' )
plt.plot(np.log(scale_factor[counter:]),np.log(k_late_radiation), label = '$n=-4$', linestyle = 'dashed' )
plt.plot(np.log(scale_factor[counter:]),np.log(k_late_f), label = '$n=-5$', linestyle = 'dashed' )
plt.plot(np.log(scale_factor[counter:]),np.log(k_late_scalar), label = '$n=-6$', linestyle = 'dashed' )
plt.ylabel('$K(a(t))$')
plt.xlabel('$a(t)$')
plt.legend()
plt.savefig('plots/Kinetic_2_log.png')

# Graph Out
plt.figure(16)
plt.plot(time_data[:counter],Viol_array[:counter])
plt.title('$\\mathcal{H}$ vs time')
plt.ylabel('$\\mathcal{H}$')
plt.xlabel('Relax Time')
plt.savefig('plots/HRelRelax.png')

plt.figure(17)
plt.plot(time_data[counter:],phi_center_data[counter:])
plt.title('$\\phi_{center}$ vs time')
plt.ylabel('$\\phi$')
plt.xlabel('Time')
plt.savefig('plots/phic.png')

print 'Program Time: ' ,time.time() - program_time, ' Seconds'
