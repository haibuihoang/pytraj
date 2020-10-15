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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.ndimage import gaussian_filter1d
from matplotlib import cm

trabw_prefix="traj_bw"
trafw_prefix="traj_fw"
out_prefix="traj_wcb"  #both backward and forward
total_traj=28
output_dir="./outputs"
data_dir="/Data/gfi/spengler/haibui/WRF_RES/ctl_lagranto"
cyclone_track="ctl_lagranto.cyclone"
dXY = 20.;

#1. read track data
df = pd.read_csv(cyclone_track, sep="\s+")
jcs=df['jsmo']
ics=df['ismo']
track_starttime=datetime.datetime.strptime( df['Time'][0], '%Y-%m-%d %H:%M:%S') 
track_2ndtime=datetime.datetime.strptime( df['Time'][1], '%Y-%m-%d %H:%M:%S') 
track_dsec=(track_2ndtime-track_starttime).total_seconds()  #Deta second

#Read backward traj and sort by time
all_traj=[]
#for i in range(1):
for i in range(total_traj):
   bwf = output_dir+"/%s.%d"%(trabw_prefix,i)
   fwf = output_dir+"/%s.%d"%(trafw_prefix,i)
   dfbw = pd.read_csv(bwf, sep="\s+")[::-1].reset_index(drop=True)
   dffw = pd.read_csv(fwf, sep="\s+")
   # combine trajectories
   df=pd.concat([dfbw,dffw[1:]],ignore_index=True)
   all_traj.append(df)




#interpolate cyclone track postion
traj=all_traj[0]
traj_times=traj['Time']
xc_interp=[]
yc_interp=[]
for i in range(len(traj_times)):
   curr_datetime=datetime.datetime.strptime(traj_times[i], '%Y-%m-%d %H:%M:%S')
   curr_secs = (curr_datetime-track_starttime).total_seconds()
   it1=math.floor(curr_secs/track_dsec)
   it2=it1+1
   t1=it1*track_dsec
   t2=it2*track_dsec
   ic=interp_1d(ics[it1],ics[it2],t1,t2,curr_secs)
   jc=interp_1d(jcs[it1],jcs[it2],t1,t2,curr_secs)
   xc_interp.append(ic*dXY)
   yc_interp.append(jc*dXY)


#Cyclone relative trajectories
for i in range(len(all_traj)):
   all_traj[i]['rx']=all_traj[i]['x']-xc_interp
   all_traj[i]['ry']=all_traj[i]['y']-yc_interp


#output the trajectories in wcb
iout=0
for ip in range(len(all_traj)):
   traj = all_traj[ip]
   dP = np.max(traj['P']) - np.min(traj['P'])
   print("Traj #%d, dP=%d"%(i,dP))
   if (dP>0):
     df=pd.DataFrame(all_traj[ip])
     outfile = output_dir+"/%s.%d"%(out_prefix,iout)
     df.to_csv(outfile,sep=' ',na_rep='-',float_format="%.4f")
     iout+=1
print("Total:",iout)


