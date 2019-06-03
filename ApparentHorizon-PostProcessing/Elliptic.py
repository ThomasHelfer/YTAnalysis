import matplotlib
#matplotlib.use('Agg')
import numpy as np
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import numpy as np


def readin(myfile):
	data = []
	myfile = open(myfile, 'r') 
	for line in myfile:
		try:
			val = float(line)	
			data.append(val)
		except: 
			donthing=0;
	return np.array(data)

N = 48
size = 0.4

minF = []
maxF = []
Cycle = []
for i in range(1,N): 
	FF = readin('Fvals%03d.csv'%i)
	uu = readin('uuvals%03d.csv'%i)
	vv = readin('vvvals%03d.csv'%i)

	dv = (np.sort(np.unique(vv))[1])	
	du = (np.sort(np.unique(uu))[1])	

	fig = plt.figure(figsize=(10,10))

	ax = fig.add_subplot(111, projection='3d')

	x = FF**2*np.sin(uu)*np.cos(vv);
	y = FF**2.*np.sin(uu)*np.sin(vv);
	z = FF**2.*np.cos(uu);
	
	Cycle.append(i)
	minF.append(min(FF**2))
	maxF.append(max(FF**2))

	# Plot the surface
	ax.scatter(x,y,z,color = 'black')

	ax.set_xlim([-size,size])
	ax.set_ylim([-size,size])
	ax.set_zlim([-size,size])

	plt.savefig("data%03d.png"%i,bbox_inches = 'tight', dpi=150)
	plt.close()

fig = plt.figure(figsize=(10,10))
plt.plot(Cycle,maxF, label = 'maxF')
plt.plot(Cycle,minF, label = 'minF')
plt.legend()
plt.savefig("Radius.png",bbox_inches = 'tight')
plt.close()

print(Cycle)
#print(maxF*np.sqrt(du*dv))
print(maxF)

fig = plt.figure(figsize=(10,10))
plt.plot(Cycle,np.array(maxF)*np.sqrt(du*dv), label = 'Estimated dx ')
plt.legend()
plt.savefig("Resolution.png",bbox_inches = 'tight')
plt.close()



#Take from 0806.4007 eq 48
def Elliptic(z):
	inte = 0
	dtheta = 0.00001
	N = int(np.floor((np.pi/2)/dtheta))
	for i in range(N):
		theta = i*dtheta
		inte += np.sqrt(1-z*(np.sin(theta))**2.0)*dtheta	
	return inte;

# a = dimensionless spin

def ratio(a):
	a = np.where(a>1,1,a)
	rplus = 1 + np.sqrt(1-a**2.0)
	ratio = np.sqrt(2*rplus)/np.pi*Elliptic((a**2/(2*rplus)))
	print("Evaluating for spin ... ") 
	return ratio 


area = np.loadtxt('area.out')
polar = np.loadtxt('polar_horizon_length.out')
equator = np.loadtxt('equator_horizon_length.out')


plt.figure(figsize=(12,8))
plt.scatter(equator[:,0],equator[:,1],label = "equator ", facecolor='red')
plt.ylabel("equator ")
plt.xlabel(r't $[1/m]$')
plt.legend(loc= "best")
plt.savefig("equator.png",bbox_inches = 'tight')
plt.close()

plt.figure(figsize=(12,8))
plt.scatter(area[:,0],area[:,1],label = "area ", facecolor='red')
plt.ylabel("area ")
plt.xlabel(r't $[1/m]$')
plt.legend(loc= "best")
#plt.ylim([34,30])
plt.savefig("area.png",bbox_inches = 'tight')
plt.close()

ratio_data = polar
ratio_data[:,1] = ratio_data[:,1]/equator[:,1]
guess = 0.74

time = area[:,0]
spin_data_1 = []
spin_data_2 = []
for i in range(len(ratio_data[:,0])):
	def search(a):
		return ratio(a)-ratio_data[i,1]
	sol = optimize.root(search,guess)
	spin_1 = sol.x
	spin_data_1.append(spin_1[0])
	guess = spin_1[0]
	print("spin1",guess)
	spin_2 = np.sqrt(1.0-((2.*np.pi*area[i,1]/(equator[i,1]**2.))-1.)**2.)
	spin_data_2.append(spin_2)
	print('--------------------------------------------------')

plt.figure(figsize=(12,8))
plt.scatter(time,spin_data_1,label = "Method 1 ", facecolor='black')
plt.scatter(time,spin_data_2,label = "Method 2 ", facecolor='red')
plt.ylabel("dimensionless spin a ")
plt.xlabel(r't $[1/m]$')
plt.legend(loc= "best")
plt.savefig("Refs.png",bbox_inches = 'tight')
plt.close()


