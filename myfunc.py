import numpy as np
import math
# Time line interpolation
# f1 and f2 are two numpy array
def interp_1d(f1,f2,t1,t2,t):
  f = f1 + (f2-f1)*(t-t1)/(t2-t1)
  return np.array(f)

#x,y,z are 1d array with the same dimension
#Var,Z has dimension of [Nz,Ny,Nx]
def interp_3d(Dx,Z,Var,x,y,z):
  num=len(x)
  pVar=np.empty_like(x)
  for ip in range(num):
    i=x[ip]/Dx
    j=y[ip]/Dx
    i1=math.floor(i)
    i2=i1+1
    j1=math.floor(j)
    j2=j1+1
    Varyz = Var[:,:,i1] + (Var[:,:,i2]-Var[:,:,i1])*(i-i1)
    Varz = Varyz[:,j1] + (Varyz[:,j2]-Varyz[:,j1])*(j-j1)
    Zyz = Z[:,:,i1] + (Z[:,:,i2]-Z[:,:,i1])*(i-i1)
    Zz = Zyz[:,j1] + (Zyz[:,j2]-Zyz[:,j1])*(j-j1)
    #Horizontal interpolation
    #Varyz=interp_1d(Var[:,:,i1],Var[:,:,i2],i1,i2,i)
    #Varz =interp_1d(Varyz[:,j1],Varyz[:,j2],j1,j2,j)
    #Zyz=interp_1d(Z[:,:,i1],Z[:,:,i2],i1,i2,i)
    #Zz =interp_1d(Zyz[:,j1],Zyz[:,j2],j1,j2,j)
    pVar[ip]=np.interp(z[ip],Zz,Varz)

  return pVar

#Var has dimension of [Nz,Ny,Nx]
def interp_2d(Dx,Var,x,y):
  num=len(x)
  pVar=np.empty_like(x)
  for ip in range(num):
    i=x[ip]/Dx
    j=y[ip]/Dx
    i1=math.floor(i)
    i2=i1+1
    j1=math.floor(j)
    j2=j1+1
    Varyz = Var[:,i1] + (Var[:,i2]-Var[:,i1])*(i-i1)
    pVar[ip] = Varyz[j1] + (Varyz[j2]-Varyz[j1])*(j-j1)
  return pVar



