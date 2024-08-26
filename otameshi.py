# Add 3 lines by S. Wakita (07Jun2017)
import matplotlib as mpl
#mpl.use('Cairo') #for png, ps, pdf, svg, ...
mpl.use('Agg') #for png
import numpy as np
import matplotlib.pyplot as plt
import pySALEPlot as psp
from sympy import *
import math
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
import os
import csv

# Open data file
model=psp.opendatfile('./Ice/jdata.dat')

# Set the distance units to cm
model.setScale('mm')

# Read the experimental data
#rtime, rexpt = np.genfromtxt('./Plotting/experimental.data', usecols=(0,1), unpack=True)
#dtime, dexpt = np.genfromtxt('./Plotting/experimental.data', usecols=(2,3), unpack=True, skip_footer=15)

Grid_size = model.inputDict['GRIDSPC']
Grid_H = model.inputDict['GRIDH'][1]
Grid_V = model.inputDict['GRIDV'][1]
CPPR    = model.inputDict['OBJRESH']
vimp    = model.inputDict['OBJVEL']
Rp      = Grid_size*CPPR
vp      = -vimp
num_impactor = model.tru[0].end
num_total = model.tracer_num
ts      = 2.*Rp/vp
dt      = -model.inputDict['DTSAVE']*ts
TEND    = -model.inputDict['TEND']
DTSAVE  = -model.inputDict['DTSAVE']
STEP_num = TEND/DTSAVE
ts_step  = 1/DTSAVE

fig = plt.figure(figsize=(12,10))
ax4 = fig.add_subplot(111)

#read No. 1
step0   = model.readStep(['TrP','TrT'],0)
step1   = model.readStep(['TrP','TrT'],model.laststep)

x1 = np.where(step0.xmark <= Grid_size*1e3/2)
x2 = x1[0]
x2_1 = x2[x2  > num_impactor ]
x3 = sorted(x2_1,reverse=True)

x3x = step0.xmark[x3]
x3y = -step0.ymark[x3]
x3P = step1.TrP[x3]*1e-9
maxp1 = max(x3P)
maxP1 =round(maxp1,1)
d_r=x3y/(Rp*1e3)

#read No. 2
model2=psp.opendatfile('../KatoIce14-2/Ice/jdata.dat')
model2.setScale('mm')
CPPR2    = model2.inputDict['OBJRESH']
Grid_size2 = model2.inputDict['GRIDSPC']
Rp2      = Grid_size2*CPPR2
num_impactor_2 = model2.tru[0].end
num_total_2 = model2.tracer_num

step0_2   = model2.readStep(['TrP','TrT'],0)
step1_2   = model2.readStep(['TrP','TrT'],model2.laststep)

x1_2 = np.where(step0_2.xmark <= Grid_size2*1e3/2)
x2_2 = x1_2[0]
x2_1_2 = x2_2[x2_2  > num_impactor_2 ]
x3_2 = sorted(x2_1_2,reverse=True)

x3x_2 = step0_2.xmark[x3_2]
x3y_2 = -step0_2.ymark[x3_2]
x3P_2 = step1_2.TrP[x3_2]*1e-9
maxp1_2 = max(x3P_2)
maxP1_2 =round(maxp1_2,1)
d_r_2=x3y_2/(Rp2*1e3)

#read No. 3
model3=psp.opendatfile('../KatoIce14-3/Ice/jdata.dat')
model3.setScale('mm')
CPPR3    = model3.inputDict['OBJRESH']
Grid_size3 = model3.inputDict['GRIDSPC']
Rp3      = Grid_size3*CPPR3
num_impactor_3 = model3.tru[0].end
num_total_3 = model3.tracer_num

step0_3   = model3.readStep(['TrP','TrT'],0)
step1_3   = model3.readStep(['TrP','TrT'],model3.laststep)

x1_3 = np.where(step0_3.xmark <= Grid_size3*1e3/2)
x2_3 = x1_3[0]
x2_1_3 = x2_3[x2_3  > num_impactor_3 ]
x3_3 = sorted(x2_1_3,reverse=True)

x3x_3 = step0_3.xmark[x3_3]
x3y_3 = -step0_3.ymark[x3_3]
x3P_3 = step1_3.TrP[x3_3]*1e-9
maxp1_3 = max(x3P_3)
maxP1_3 =round(maxp1_3,1)
d_r_3=x3y_3/(Rp3*1e3)



x4 = np.where(step0.xmark <= Grid_size/2*1e3)
x5 = x4[0]
x6 = x5[x5  < num_impactor ]
x6x = step0.xmark[x6]
x6y = step0.ymark[x6]
x6P = step1.TrP[x6]*1e-9
x6T = step1.TrT[x6]
maxp61 = max(x6P[0:20])
maxt61 = max(x6T[0:20])
maxP61 =round(maxp61,1)
maxT61 =round(maxt61,0)
maxp62 = max(x6P)
maxt62 = max(x6T)
maxP62 =round(maxp62,1)
maxT62 =round(maxt62,0)

prelabelscale = len(str(maxp1).split('.')[0])
prelabel1 = 10**(prelabelscale)
prelabel2 = prelabel1/2

if maxp1 <= prelabel2:
        prelabel1 = prelabel2
        prelabel2 = float(prelabel2)/2

maxy1 = max(d_r)

pre = []
dt_r   = []
impact_point = 0
for i in np.arange(0,model.laststep-1,1):
#for i in np.arange(0,20,1):
  step   = model.readStep(['Trp','TrP'],i)
  dt_t = dt*i
  pre.append(step.Trp[0]*1e-9)
  dt_r.append(dt_t)

#cb1.remove()
#ax2.remove()
# Plot the pressure field
for u in range(model.tracer_numu): #loop over tracer clouds
         tstart = model.tru[u].start
         tend = model.tru[u].end
         scat1 = ax4.scatter(step0.xmark[tstart:tend],step0.ymark[tstart:tend],
                c=step1.TrP[tstart:tend]*1e-6,vmin=38,vmax=1e+3,
                 cmap='viridis',s=40,linewidths=0)

cb1=fig.colorbar(scat1,ax=ax4)
cb1.set_label('Tracer Prak Pressure [MPa]')
ax4.set_xlabel('x  [mm]')
ax4.set_ylabel('z  [mm]')
ax4.set_xlim([0,8*Rp*1e3])
ax4.set_ylim([-10*Rp*1e3,3*Rp*1e3])


x4 = np.where(Grid_size*1e3/2 <= step0.xmark)
x5 = x4[0]
x6 = x5[num_impactor <= x5 ]
x6x = step0.xmark[x6]
x6y = -step0.ymark[x6]
x6P = step1.TrP[x6]
maxp3 = max(x6P)

# x_g = len(x6)/(Grid_H-1)
y_g = int(len(x6)/(Grid_V-1))
x_g = len(x6)/(Grid_H)

y_range = (x_g-CPPR*10)
x_range = (CPPR*10)
x10 =[]

for i in np.arange(0,x_range,1):
  x7 = x6[step0.xmark[x6] <= (i+1)*Grid_size*1e3]
  x8 = x7[i*Grid_size*1e3<= step0.xmark[x7] ]
  x9 = x8[y_range:]
  x9P = step1.TrP[x9]

  for j in np.arange(0,len(x9)-1,1):
    if x9P[j] <= 38e+6 and x9P[j+1] >= 38e+6:
      thre_id = num_impactor+i*x_g+j+y_range
      x10.append(thre_id)
      print('bottom',thre_id,i,j,step0.xmark[thre_id],step0.ymark[thre_id])
    if x9P[j] >= 38e+6 and x9P[j+1] <= 38e+6:
      thre_id = num_impactor+i*x_g+j+y_range
      x10.append(thre_id)
      print('top',thre_id,i,j,step0.xmark[thre_id],step0.ymark[thre_id])

ax4.scatter(step0.xmark[x10],step0.ymark[x10],c='white',s=5)

#fig.tight_layout(rect=[0,0,0.98,0.9])
# Save the figure
cwd = os.getcwd()
dirname = os.path.basename(cwd)

plt.savefig('../property_iso/{}isobaric_2.png'.format(dirname))
