# purposes:
# Combine backward and forward trajectories
# filter out warm conveyor belt
# add cyclone-relative coordinates
import numpy as np
import netCDF4
from wrf import getvar
import datetime
import math
from myfunc import *
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib import cm
#import fix_axes3d

tra_prefix="traj_wcb.400m"  #both backward and forward
figout="wcb.400m.pdf"
total_traj=51
output_dir="./outputs"
data_dir="/Data/gfi/spengler/haibui/WRF_RES/ctl_lagranto"
dXY = 20.;
base_date=datetime.datetime(1001,1,1)

#for i in range(1):
all_traj=[]
for i in range(total_traj):
   f = output_dir+"/%s.%d"%(tra_prefix,i)
   df = pd.read_csv(f, sep="\s+")
   print(df.min()['P'])
   # filter here
   if ((df.max()['P']-df.min()['P'])>55000.):
     all_traj.append(df)
     n_points=len(df)
it0 = 990
num_traj=len(all_traj)

traj_start=datetime.datetime.strptime(all_traj[0]['Time'][0], '%Y-%m-%d %H:%M:%S')
deltat=datetime.datetime.strptime(all_traj[0]['Time'][1], '%Y-%m-%d %H:%M:%S')-traj_start
hh0=(traj_start-base_date).seconds/3600.
dsec=deltat.seconds
print(hh0,dsec)
hh=np.empty([len(all_traj[0]['Time'])],dtype=float)
for i in range(len(all_traj[0]['Time'])):
  hh[i]= hh0 + i*dsec/3600.

print(hh)


traj=all_traj[0]
curr_datetime=datetime.datetime.strptime(traj['Time'][it0], '%Y-%m-%d %H:%M:%S')
print("Reading ",curr_datetime)
wrf = netCDF4.Dataset(data_dir+"/wrfout_d01_"+curr_datetime.strftime("%Y-%m-%d_%H:%M:%S"))
xc0=traj['x'][it0]-traj['rx'][it0]
yc0=traj['y'][it0]-traj['ry'][it0]
ic0=int(xc0/dXY)
jc0=int(yc0/dXY)
Hdom = 1000. #Half domain to be plot
dij = int(Hdom/dXY)
(i1,i2,j1,j2)=(ic0-dij,ic0+dij,jc0-dij,jc0+dij)

slp =  np.array(getvar(wrf,"slp"))
slpm = slp.mean(axis=1)
slpa = np.subtract(slp,slpm.reshape((len(slpm),1)))[j1:j2+1,i1:i2+1] 


xarr = np.arange(-Hdom,Hdom+dXY,dXY,dtype=float)
yarr = xarr
xx,yy = np.meshgrid(xarr, yarr)


#-----------------------
#Plotting....
print("Ploting....")
#cyclone center at t0




fig = plt.figure(figsize=(10, 10))
gs = mpl.gridspec.GridSpec(5,1)

projection='3d'
ax = fig.add_subplot(gs[0:3,0], projection='3d')
ax1 = fig.add_subplot(gs[3,0])
ax2 = fig.add_subplot(gs[4,0])

rb = plt.get_cmap('gist_rainbow')

for i in range(len(all_traj)):
  color=rb(1.*i/num_traj)
  traj=all_traj[i]
  ax1.plot(hh,traj['QV'],c=color,linewidth=0.5)

for i in range(len(all_traj)):
  color=rb(1.*i/num_traj)
  traj=all_traj[i]
  ax2.plot(hh,traj['P'],c=color,linewidth=0.5)



#ax = Axes3D(fig,proj_type="ortho")
#ax.get_proj = lambda: 1.5*np.dot(Axes3D.get_proj(ax), np.diag([1, 0.5, 0.5, 1]))
ax.grid(False)
ax.set_zlim(0., 7.)
ax.set_xlim(-Hdom, 2*Hdom)
ax.set_ylim(-Hdom, Hdom)
ax.view_init(elev=15., azim=-80)

ax.contourf(xx, yy, slpa, zdir='z', offset=0, cmap=cm.Blues_r,zorder=-99)
ax.contour(xx, yy, slpa, colors="black",linewidths=0.5,zdir='z', offset=0,zorder=-90)
fig.tight_layout()

for i in range(len(all_traj)):
   traj=all_traj[i]
   dP = np.max(traj['P']) - np.min(traj['P'])
   segments=[]
   for j in range(len(traj['rx'])-1):
       segments.append([ [traj['rx'][j],traj['ry'][j],traj['z'][j]],[traj['rx'][j+1],traj['ry'][j+1],traj['z'][j+1]] ])
   segments=np.array(segments)
   lc = Line3DCollection(segments, cmap=plt.get_cmap('turbo'),
                    norm=plt.Normalize(8, 14))
   lc.set_array(traj['QV'][1:])
   lc.set_linewidth(2)
   lc.set_zorder(300)
   line=ax.plot(traj['rx'], traj['ry'], 0.1,c="gray",linewidth=1,alpha=0.5,zorder=20)
  
   ax.add_collection3d(lc)


fig.colorbar(lc, ax=ax)
plt.savefig(figout)

plt.show()

