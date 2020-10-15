#perform backward trajectory analysis of WRF-ARW idealized output
import numpy as np
import netCDF4 
from wrf import getvar
import datetime
import math
from myfunc import *
import pandas as pd
import time
import matplotlib.pyplot as plt

start_z=600  # in m, positive value will overide z in start.xyz

data_dir="/Data/gfi/spengler/haibui/WRF_RES/ctl_lagranto"
output_dir="./outputs"
base_date=datetime.datetime(1001,1,1)
output_freq=1 #How many hours between file outputs
start_hours=24
end_hours=72
output_prefix="traj_fw"  #fw,bw
dt = 120  # in second, end_hours<start_hoursmeans backward integration
max_iteration=5 #Iteration of lagrangian integration 
if (start_z>0):
  output_prefix=output_prefix+".%dm"%start_z
  print(output_prefix)

#physical parameters
fCor=1.028e-4
Rd = 287
RdCp = 0.286
P0=100000.
dXY=20000.

start_date=base_date+datetime.timedelta(hours=start_hours)
stop_date=base_date+datetime.timedelta(hours=end_hours)


print("Reading data at "+start_date.strftime("%Y-%m-%d_%H:%M:%S"))
sfile = netCDF4.Dataset(data_dir+"/wrfout_d01_"+start_date.strftime("%Y-%m-%d_%H:%M:%S"))

Ut = np.array(getvar(sfile,"ua"))
Vt = np.array(getvar(sfile,"va"))
Wt = np.array(getvar(sfile,"wa"))
Zt = np.array(getvar(sfile,"z"))
Pt = np.array(getvar(sfile,"pres"))
RHt = np.array(getvar(sfile,"rh"))
TKt= np.array(getvar(sfile,"tk"))
QVt = sfile['QVAPOR'][0,:,:,:]  #Convert to g/Kg
Rhot = Pt/(Rd*TKt*(1+0.61*QVt))
THt = TKt*np.power((P0/Pt),RdCp)
Vor = np.gradient(Vt,dXY,axis=2) - np.gradient(Ut,dXY,axis=1)
dZ = np.gradient(Zt,axis=0)
dTH = np.gradient(THt,axis=0)
PVt = 1e+6*( fCor + Vor )*dTH/dZ/Rhot #In;PVU
LHFt = sfile['LH'][0,:,:]  #W.m-2
SHFt = sfile['HFX'][0,:,:]  #W.m-2


[Nz,Ny,Nx]=Ut.shape
Dx=sfile.DX
print("--------------------------")
print("  WRF information")
print("     Nx: ",Nx)
print("     Ny: ",Ny)
print("     Nz: ",Nz)
print("     Dx: ",Dx)
Domaininfo={"Nx":Nz,"Ny":Ny,"Nz":Nz,"Dx":Dx}

curr_time=start_hours*3600.
end_time=end_hours*3600.

#read start position
df = pd.read_csv('start.xyz', sep="\s+")
print("--------------------------")
print("  Start time: "+start_date.strftime("%Y-%m-%d_%H:%M:%S"))
print("  Stop time : "+stop_date.strftime("%Y-%m-%d_%H:%M:%S"))

if (end_hours<start_hours):
  dt=-dt  # Negative time step for backward integration
  print("  Backward integration")
  (U1,V1,W1,Z1,P1,RH1,QV1,TH1,PV1,LHF1,SHF1) = (Ut,Vt,Wt,Zt,Pt,RHt,QVt,THt,PVt,LHFt,SHFt)
  hh1_old=start_hours
else:
  print("  Forward integration")
  (U2,V2,W2,Z2,P2,RH2,QV2,TH2,PV2,LHF2,SHF2)  = (Ut,Vt,Wt,Zt,Pt,RHt,QVt,THt,PVt,LHFt,SHFt)
  hh1_old=start_hours-output_freq


if (start_z>0.):
  df['z'] = start_z/1000.

print("  Lagrangian timestep:",dt)
print("  Start positions: ")
num_particles=df.shape[0]

print(df)
#Convert from km to m
px=np.array(df['x'],dtype=np.float)*1000.
py=np.array(df['y'],dtype=np.float)*1000.
pz=np.array(df['z'],dtype=np.float)*1000.


pu=interp_3d(Dx,Zt,Ut,px,py,pz)
pv=interp_3d(Dx,Zt,Vt,px,py,pz)
pw=interp_3d(Dx,Zt,Wt,px,py,pz)
pP=interp_3d(Dx,Zt,Pt,px,py,pz)
pTH=interp_3d(Dx,Zt,THt,px,py,pz)
pPV=interp_3d(Dx,Zt,PVt,px,py,pz)
pRH=interp_3d(Dx,Zt,RHt,px,py,pz)
pQV=interp_3d(Dx,Zt,QVt,px,py,pz)
pLHF=interp_2d(Dx,LHFt,px,py)
pSHF=interp_2d(Dx,SHFt,px,py)

# Initialize trajectories
Trajx=[px]
Trajy=[py]
Trajz=[pz]
Trajt=[start_date]
TrajP=[pP]
TrajPV=[pPV]
TrajRH=[pRH]
TrajQV=[pQV]
TrajTH=[pTH]
TrajLHF=[pLHF]
TrajSHF=[pSHF]

Integration=True
while (Integration):
   #Define condtion for the integration depend on backward or forward
   curr_time = curr_time+dt   # This actually the new time step
   if (dt<0):
      Integration=(curr_time > end_time)
   else:
      Integration=(curr_time < end_time)

   curr_date = base_date+datetime.timedelta(seconds=curr_time)
   hh1 = math.floor(curr_time/(3600*output_freq))
   hh2 = hh1+output_freq

   print("Current time: %4.2fh"%(curr_time/3600.)),
     
   # Take quite long time to read data!!!
   # We should read data only when neccessary!!
   # And we should only read one time step
   if (hh1 != hh1_old):
     if (dt<0):
        date1 = base_date+datetime.timedelta(hours=hh1)
        read_date=date1
     else:
        date2 = base_date+datetime.timedelta(hours=hh2)
        read_date=date2

     print("Reading data at "+start_date.strftime("%Y-%m-%d_%H:%M:%S"))
     sfile = netCDF4.Dataset(data_dir+"/wrfout_d01_"+read_date.strftime("%Y-%m-%d_%H:%M:%S"))
     Ut = np.array(getvar(sfile,"ua"))
     Vt = np.array(getvar(sfile,"va"))
     Wt = np.array(getvar(sfile,"wa"))
     Zt = np.array(getvar(sfile,"z"))
     Pt = np.array(getvar(sfile,"pres"))
     RHt = np.array(getvar(sfile,"rh"))
     TKt= np.array(getvar(sfile,"tk"))
     QVt = sfile['QVAPOR'][0,:,:,:]
     Rhot = Pt/(Rd*TKt*(1+0.61*QVt))
     THt = TKt*np.power((P0/Pt),RdCp)
     Vor = np.gradient(Vt,dXY,axis=2) - np.gradient(Ut,dXY,axis=1)
     dZ = np.gradient(Zt,axis=0)
     dTH = np.gradient(THt,axis=0)
     PVt = 1e+6*( fCor + Vor )*dTH/dZ/Rhot #In;PVU
     LHFt = sfile['LH'][0,:,:]  #W.m-2
     SHFt = sfile['HFX'][0,:,:]  #W.m-2


     if (dt<0):
      (U2,V2,W2,Z2,P2,RH2,QV2,TH2,PV2,LHF2,SHF2) = (U1,V1,W1,Z1,P1,RH1,QV1,TH1,PV1,LHF1,SHF1)
      (U1,V1,W1,Z1,P1,RH1,QV1,TH1,PV1,LHF1,SHF1) = (Ut,Vt,Wt,Zt,Pt,RHt,QVt,THt,PVt,LHFt,SHFt)
      #read from file 1
     else:
      (U1,V1,W1,Z1,P1,RH1,QV1,TH1,PV1,LHF1,SHF1) = (U2,V2,W2,Z2,P2,RH2,QV2,TH2,PV2,LHF2,SHF2)
      (U2,V2,W2,Z2,P2,RH2,QV2,TH2,PV2,LHF2,SHF2)  = (Ut,Vt,Wt,Zt,Pt,RHt,QVt,THt,PVt,LHFt,SHFt)

   U = interp_1d(U1,U2,hh1*3600.,hh2*3600.,curr_time) 
   V = interp_1d(V1,V2,hh1*3600.,hh2*3600.,curr_time) 
   W = interp_1d(W1,W2,hh1*3600.,hh2*3600.,curr_time) 
   Z = interp_1d(Z1,Z2,hh1*3600.,hh2*3600.,curr_time) 
   P = interp_1d(P1,P2,hh1*3600.,hh2*3600.,curr_time) 
   PV = interp_1d(PV1,PV2,hh1*3600.,hh2*3600.,curr_time) 
   TH = interp_1d(TH1,TH2,hh1*3600.,hh2*3600.,curr_time) 
   RH = interp_1d(RH1,RH2,hh1*3600.,hh2*3600.,curr_time) 
   QV = interp_1d(QV1,QV2,hh1*3600.,hh2*3600.,curr_time) 
   LHF = interp_1d(LHF1,LHF2,hh1*3600.,hh2*3600.,curr_time) 
   SHF = interp_1d(SHF1,SHF2,hh1*3600.,hh2*3600.,curr_time) 
   hh1_old=hh1
   hh2_old=hh2

   for i in range(max_iteration):
      if (i==0):
        pxnew = px + pu*dt
        pynew = py + pv*dt
        pznew = pz + pw*dt
      else:
        pxnew = px + 0.5*(pu+punew)*dt
        pynew = py + 0.5*(pv+pvnew)*dt
        pznew = pz + 0.5*(pw+pwnew)*dt

      pznew = np.where(pznew < 0, 0, pznew)  # Limit the boundary of z
      punew=interp_3d(Dx,Z,U,pxnew,pynew,pznew)
      pvnew=interp_3d(Dx,Z,V,pxnew,pynew,pznew)
      pwnew=interp_3d(Dx,Z,W,pxnew,pynew,pznew)
   

   #Now interpolate for other fields
   TrajP.append(interp_3d(Dx,Z,P,pxnew,pynew,pznew))
   TrajPV.append(interp_3d(Dx,Z,PV,pxnew,pynew,pznew))
   TrajTH.append(interp_3d(Dx,Z,TH,pxnew,pynew,pznew))
   TrajRH.append(interp_3d(Dx,Z,RH,pxnew,pynew,pznew))
   TrajQV.append(interp_3d(Dx,Z,QV,pxnew,pynew,pznew))
   TrajLHF.append(interp_2d(Dx,LHF,pxnew,pynew))
   TrajSHF.append(interp_2d(Dx,SHF,pxnew,pynew))
     
   Trajx.append(pxnew)
   Trajy.append(pynew)
   Trajz.append(pznew)
   Trajt.append(curr_date)
  
   (pu,pv,pw)=(punew,pvnew,pwnew)
   (px,py,pz)=(pxnew,pynew,pznew)

   #endt = time.time()
   #elapsed=endt-startt
   #print("Calculating: ",elapsed)
   

Trajx=np.array(Trajx)/1000. #Convert to km
Trajy=np.array(Trajy)/1000.
Trajz=np.array(Trajz)/1000.
TrajP=np.array(TrajP)
TrajPV=np.array(TrajPV)
TrajTH=np.array(TrajTH)
TrajRH=np.array(TrajRH)
TrajQV=np.array(TrajQV)
TrajLHF=np.array(TrajLHF)
TrajSHF=np.array(TrajSHF)
#output trajectories to files
for ip in range(num_particles):
   df=pd.DataFrame({'Time':Trajt,'x':Trajx[:,ip],'y':Trajy[:,ip],'z':Trajz[:,ip],
                    'P':TrajP[:,ip],'PV':TrajPV[:,ip],'TH':TrajTH[:,ip],'RH':TrajRH[:,ip],'QV':1000.*TrajQV[:,ip],'LHF':TrajLHF[:,ip],'SHF':TrajSHF[:,ip]})
   outfile = output_dir+"/%s.%d"%(output_prefix,ip)
   df.to_csv(outfile,sep=' ',na_rep='-',float_format="%.4f")



