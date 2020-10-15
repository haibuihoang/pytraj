"""
Read and write a time SLP, Theta850 of WRF
"""
import numpy as np
import netCDF4 as nc4
from wrf import getvar,interplevel
import datetime
import matplotlib.pyplot as plt

ifile="/Data/gfi/spengler/haibui/WRF_RES/ctl_lagranto/wrfout_d01_1001-01-03_03:00:00"
ofile="wrftmp.nc"

wrf = nc4.Dataset(ifile)
slp =  np.array(getvar(wrf,"slp"))
slpm = slp.mean(axis=1)

p = getvar(wrf, "pressure")
th = getvar(wrf, "th")

th_850=interplevel(th, p, 850.)
#plt.contour(th_850)
#plt.show()


dim=slp.shape
ny=dim[0]
nx=dim[1]


f = nc4.Dataset(ofile,'w', format='NETCDF4')
#grp = f.createGroup('data')
f.createDimension('x', nx)
f.createDimension('y', ny)
slpa = f.createVariable('slpa', 'f4', ('y', 'x'))
th850 = f.createVariable('th850', 'f4', ('y', 'x'))
slpa[:,:] = np.subtract(slp,slpm.reshape((len(slpm),1)))
th850[:,:]=th_850
f.close()

