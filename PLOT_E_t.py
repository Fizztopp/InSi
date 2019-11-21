import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

tc = 0.658 # 1/eV = 0.658 fs        

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 15
mpl.rcParams['font.size'] = 20  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['axes.titlesize'] = 24
mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 28
mpl.rcParams['figure.figsize'] = [10.,5]

m = 1

RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 

TIMESTEPS = 1e4

file_BANDS = open('mu.txt','r')
mu  = 0
file_BANDS.close()

file_BANDS = open('ETOT_t.txt','r')
MAT_E = np.loadtxt(file_BANDS)-mu
file_BANDS.close()

file_Driving = open('DRIVING_t.txt','r')
Driving = np.loadtxt(file_Driving)
file_Driving.close()

def gauss(sigma, shift, x):
    return np.exp(-0.5*((x-shift)/sigma)**2)

facx = 188./3.84*0.1 # E_max = OMEGA*A_max = 200meV/4.20AA

t = np.linspace(0,3039.51,TIMESTEPS)*tc
fig2 = plt.figure(2)
gs2 = gridspec.GridSpec(1, 1)

ax11 = fig2.add_subplot(gs2[0,0])
ax11.set_xlabel(r'$\mathrm{time}$ $\mathrm{(fs)}$',fontsize=24)
ax11.set_ylabel(r'$\mathrm{\Delta E_{tot}}$ $\mathrm{(eV)}$',fontsize=24)
MAT_E=MAT_E[:]-MAT_E[0]   
ax11.plot(t,MAT_E, linewidth=2.0)  
ax11.set_xlim(500, 1500)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) 
plt.show()

