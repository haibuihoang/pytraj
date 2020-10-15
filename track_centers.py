# purposes:
# Given the first guesss
# track the cyclone center as the minimum of SLP-anomaly
#
import numpy as np
import numpy.ma as ma
import netCDF4
from wrf import getvar,smooth2d
import datetime
import math
from myfunc import *
import pandas as pd
import time
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

data_dir="/Data/gfi/spengler/haibui/WRF_RES/ctl_lagranto"
data_freq=1 #How many hours between file outputs
outfile="./ctl_lagranto.cyclone"

jc0,ic0 = 150,100  # first guess center
dij=20  #search 40x40 gridpoint (800x800km)

base_date=datetime.datetime(1001,1,1)
start_hours=0
end_hours=120

start_date=base_date+datetime.timedelta(hours=start_hours)
end_date=base_date+datetime.timedelta(hours=end_hours)

curr_date=start_date
ics=[]
jcs=[]
times=[]
while(curr_date<end_date):
    #do the searching
    print(curr_date)
    wrffile = netCDF4.Dataset(data_dir+"/wrfout_d01_"+curr_date.strftime("%Y-%m-%d_%H:%M:%S"))
    slp =  np.array(getvar(wrffile,"slp"))
    slpm = slp.mean(axis=1)
    slpa = np.subtract(slp,slpm.reshape((len(slpm),1))) #Need to reshape to be correctly broadcasted

    #search within
    slpaa=slpa[jc0-dij:jc0+dij,ic0-dij:ic0+dij]
    slpaa=smooth2d(slpaa,passes=5)
    (jca,ica) = np.unravel_index(slpaa.argmin(),slpaa.shape)
    (jc,ic)=(jca+jc0-dij,ica+ic0-dij)
    print((jc,ic))
    ics.append(ic)
    jcs.append(jc)
    times.append(curr_date)
    #plt.contour(slpa)
    #plt.show()
    (jc0,ic0)=(jc,ic)
    curr_date = curr_date+datetime.timedelta(hours=data_freq)

ics=np.array(ics)	
jcs=np.array(jcs)	

ismo = gaussian_filter1d(ics.astype(float), 3,mode="nearest")           
jsmo = gaussian_filter1d(jcs.astype(float), 3,mode="nearest")    

df=pd.DataFrame({'Time':times,'ic':ics,'jc':jcs,'ismo':ismo,'jsmo':jsmo})
df.to_csv(outfile,sep=' ',na_rep='-',float_format="%.4f")

