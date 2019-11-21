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
mpl.rcParams['legend.fontsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 28
mpl.rcParams['figure.figsize'] = [5.,6]

m = 1

RED = '#e41a1c'
BLUE = '#6baed6'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 


file_BANDS = open('mu.txt','r')
mu  = np.loadtxt(file_BANDS)
file_BANDS.close()


file_BANDS = open('bands0.txt','r')
MAT_BANDS = np.loadtxt(file_BANDS)+0.046433
file_BANDS.close()

N_BAND = np.size(MAT_BANDS[:,0])
k=np.linspace(-np.pi/2,np.pi/2,N_BAND)     
print(N_BAND)    

fig1 = plt.figure(2)
gs1 = gridspec.GridSpec(1, 1)
ax11 = fig1.add_subplot(gs1[0,0])
ax11.patch.set_facecolor(BLUE)
ax22 = ax11.twiny()
ax11.set_xticklabels([r'$\mathrm{X_{4 \times 2}}$', r'$\Gamma_{4 \times 2}$', r'$\mathrm{X_{4 \times 2}}$'], fontsize=20)
ax11.set_xticks([-np.pi/2, 0, np.pi/2])
ax22.set_xticklabels([r'', r'$\mathrm{X_{4 \times 1}}$', r''], fontsize=20)
ax22.xaxis.tick_top()
ax22.set_xticks([-np.pi/2, 0, np.pi/2])
ax11.set_ylabel(r'energy (eV)',fontsize=20)
#ax11.set_xlabel(r'$\mathrm{k}$',fontsize=24)
ax11.plot(k,[0]*N_BAND,'k--', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,0], 'k', linewidth=2.0, label=r"distortion")
ax11.plot(k,MAT_BANDS[:,1], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,2], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,3], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,4], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,5], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,6], 'k', linewidth=2.0)
ax11.plot(k,MAT_BANDS[:,7], 'k', linewidth=2.0)
ax11.set_ylim(-1.0,1.5)
ax11.set_xlim(-np.pi/2,np.pi/2)
         
plt.tight_layout()
plt.show()
