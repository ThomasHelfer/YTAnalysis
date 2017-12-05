#########################
# Preamble
########################

import numpy as np
import pylab as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import string
from scipy import signal

# ============================
#  Doddy defintions
# ============================

F1 = 30 # Axes font size
F2 = 15 # Legend font size
line = 2.5 # Line width
alp = 1 # alpha
scale= 'linear'
colorshift=3

# Setting the font structure

rc = mpl.rcParams # Font structure is called rc now
rc['text.usetex'] = True # Tex fonts
rc['font.family'] = 'serif'
rc['font.serif'].insert(0,'cm') # Default font is computer modern for latex
rc['font.size'] = F1
rc['xtick.labelsize'] = 'small'
rc['ytick.labelsize'] = 'small'
rc['legend.fontsize'] = F2

# Set a color stucture

# These are the "Tableau 20" colors as RGB.  
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.)

names = np.array([
'Weyl4_l2_m0_data.out',
'Weyl4_l2_m1_data.out',
'Weyl4_l2_m2_data.out',
'Weyl4_l2_m1n_data.out',
'Weyl4_l2_m2n_data.out',
'Weyl4_l3_m0_data.out',
'Weyl4_l3_m1_data.out',
'Weyl4_l3_m2_data.out',
'Weyl4_l3_m3_data.out',
'Weyl4_l3_m1n_data.out',
'Weyl4_l3_m2n_data.out',
'Weyl4_l3_m3n_data.out',
'Weyl4_l4_m0_data.out',
'Weyl4_l4_m1_data.out',
'Weyl4_l4_m2_data.out',
'Weyl4_l4_m3_data.out',
'Weyl4_l4_m4_data.out',
'Weyl4_l4_m1n_data.out',
'Weyl4_l4_m2n_data.out',
'Weyl4_l4_m3n_data.out',
'Weyl4_l4_m4n_data.out',
'Weyl4_l5_m0_data.out',
'Weyl4_l5_m1_data.out',
'Weyl4_l5_m2_data.out',
'Weyl4_l5_m3_data.out',
'Weyl4_l5_m4_data.out',
'Weyl4_l5_m5_data.out',
'Weyl4_l5_m1n_data.out',
'Weyl4_l5_m2n_data.out',
'Weyl4_l5_m3n_data.out',
'Weyl4_l5_m4n_data.out',
'Weyl4_l5_m5n_data.out',
'Weyl4_l6_m0_data.out',
'Weyl4_l6_m1_data.out',
'Weyl4_l6_m2_data.out',
'Weyl4_l6_m3_data.out',
'Weyl4_l6_m4_data.out',
'Weyl4_l6_m5_data.out',
'Weyl4_l6_m6_data.out',
'Weyl4_l6_m1n_data.out',
'Weyl4_l6_m2n_data.out',
'Weyl4_l6_m3n_data.out',
'Weyl4_l6_m4n_data.out',
'Weyl4_l6_m5n_data.out',
'Weyl4_l6_m6n_data.out',
])

labes = np.array([
'l = 2 m = 0 ',
'l = 2 m = 1 ',
'l = 2 m = 2 ',
'l = 2 m = -1 ',
'l = 2 m = -2 ',
'l = 3 m = 0 ',
'l = 3 m = 1 ',
'l = 3 m = 2 ',
'l = 3 m = 3 ',
'l = 3 m = -1 ',
'l = 3 m = -2 ',
'l = 3 m = -3 ',
'l = 4 m = 0 ',
'l = 4 m = 1 ',
'l = 4 m = 2 ',
'l = 4 m = 3 ',
'l = 4 m = 4 ',
'l = 4 m = -1 ',
'l = 4 m = -2 ',
'l = 4 m = -3 ',
'l = 4 m = -4 ',
'l = 5 m = 0 ',
'l = 5 m = 1 ',
'l = 5 m = 2 ',
'l = 5 m = 3 ',
'l = 5 m = 4 ',
'l = 5 m = 5 ',
'l = 5 m = -1 ',
'l = 5 m = -2 ',
'l = 5 m = -3 ',
'l = 5 m = -4 ',
'l = 5 m = -5 ',
'l = 6 m = 0 ',
'l = 6 m = 1 ',
'l = 6 m = 2 ',
'l = 6 m = 3 ',
'l = 6 m = 4 ',
'l = 6 m = 5 ',
'l = 6 m = 6 ',
'l = 6 m = -1 ',
'l = 6 m = -2 ',
'l = 6 m = -3 ',
'l = 6 m = -4 ',
'l = 6 m = -5 ',
'l = 6 m = -6 ',
])

DATA = np.array([
"ANALYSIS",
"ANALYSIS2",
"ANALYSIS3",
"ANALYSIS4",
#"ANALYSIS5",
#"ANALYSIS6",
#"ANALYSIS7",
#"ANALYSIS8",
#"ANALYSIS9",
#"ANALYSIS10",
#"ANALYSIS11",
#"ANALYSIS12",
#"ANALYSIS13",
#"ANALYSIS14",
#"ANALYSIS15",
#"ANALYSIS16",
#"ANALYSIS17",
#"ANALYSIS18",
#"ANALYSIS19",
#"ANALYSIS20",
#"ANALYSIS21",
#"ANALYSIS22",
])

N = len(names)
NDATA = len(DATA)

# ============================
#  Prepare data
# ============================

def Prepare(str):
	with open(str, 'r') as file :
  		filedata = file.read()
	filedata = filedata.replace('+-','-')
	# Write the file out again
	with open(str, 'w') as file:
  		file.write(filedata)

for j in range (NDATA):
	# r = 60
	
	dic = (DATA[j]+"/Psi4_r60/")

	for i in range(N): 
		Prepare((dic+names[i]))
	
#	dic = (DATA[j]+"/Psi4_r75/")
#
#	for i in range(N): 
#		Prepare((dic+names[i]))

#	dic = (DATA[j]+"/Psi4_r100/")
#
#	for i in range(N): 
#		Prepare((dic+names[i]))

# ============================
#  DATA read in 
# ============================

#r = 60 
dic = "/Psi4_r60/"

Weyldata1_60 = [[[]for _ in range(NDATA)]for _ in range(N)]

time1_60 = [[]for _ in range(NDATA)]

for j in range(NDATA):
	time1_60[j].append(np.loadtxt((DATA[j]+dic+'time.out'),unpack=True))
	for i in range(N):
		data = np.loadtxt((DATA[j]+dic+names[i]),unpack=True,dtype = "complex")
		Weyldata1_60[i][j].append(data)
#		b, a = signal.butter(8, 0.125)
#		lowpass = signal.filtfilt(b, a, data)
#		Weyldata1_60[i][j].append(lowpass)

# cleaning data
Fuseddata = [[]for _ in range(N)]
Fuseddata_cleand = [[]for _ in range(N)]
Fusetime  = [[]for _ in range(N)]
time0 = min(time1_60[0][0])
for j in range(N):
	plt.figure(figsize=(15,10))
	Weyldata1_60_temp = np.array(Weyldata1_60[j][0][0])
	time_60_temp = np.array(time1_60[0][0])
	for i in range(NDATA-1):	
		Weyldata1_60_temp = np.hstack((Weyldata1_60_temp,Weyldata1_60[j][i+1][0]))
		time_60_temp = np.hstack((time_60_temp,time1_60[i+1][0]))
	Fusetime[j].append(time_60_temp)
	Fuseddata[j].append(Weyldata1_60_temp)
	b, a = signal.butter(8, 0.125)
	lowpass = signal.filtfilt(b, a, Weyldata1_60_temp)
	Fuseddata_cleand[j].append(lowpass)
	plt.plot(Fusetime[j][0]-time0-60,np.real(Fuseddata[j][0]),linewidth = line, alpha = alp, color = "darkviolet", label = "signal")
#	plt.plot(Fusetime[j][0]-time0-60,np.real(Fuseddata_cleand[j][0]),linewidth = line*1.25, alpha = alp, color = "Green", label = "with low-pass filter ")
#	plt.xlim([0,max(Fusetime[j][0]-time0-60)])
	plt.legend([(labes[j] + " rad = 60 "),(labes[j] + " rad = 60 + filter")],loc = "best")
	plt.xlabel(r'$t_{ret}~[1/m_{a}]$')
	plt.ylabel(r'$r\Psi_4$')
	plt.grid()
	plt.savefig(("Weyl_cleaned_("+labes[j]+')_real'+'.png'),bbox_inches = 'tight')
	plt.close()

















