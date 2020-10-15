from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np

Axes3D.scale_x=1
Axes3D.scale_y=1
Axes3D.scale_z=1
def get_proj_scale(self):
    """                                                                                                  Solution from: https://stackoverflow.com/questions/30223161/matplotlib-mplot3d-how-to-increase-the-size-of-an-axis-stretch-in-a-3d-plo
    """
    #Axes are scaled down to fit in scene                                                                                                                
    max_scale=max(self.scale_x, self.scale_y, self.scale_z)

    scale_x=self.scale_x/max_scale
    scale_y=self.scale_y/max_scale
    scale_z=self.scale_z/max_scale

    #Create scaling matrix                                                                                                                                                                                                                   
    scale = np.array([[scale_x,0,0,0],
                      [0,scale_y,0,0],
                      [0,0,scale_z,0],
                      [0,0,0,1]])
    #print(scale)    
    
    relev, razim = np.pi * self.elev/180, np.pi * self.azim/180

    xmin, xmax = self.get_xlim3d()
    ymin, ymax = self.get_ylim3d()
    zmin, zmax = self.get_zlim3d()

    # transform to uniform world coordinates 0-1.0,0-1.0,0-1.0                                                                                                                                                                             
    worldM = proj3d.world_transformation(
        xmin, xmax,
        ymin, ymax,
        zmin, zmax)

    # look into the middle of the new coordinates                                                                                                                                                                                          
    R = np.array([0.5, 0.5, 0.5])

    xp = R[0] + np.cos(razim) * np.cos(relev) * self.dist
    yp = R[1] + np.sin(razim) * np.cos(relev) * self.dist
    zp = R[2] + np.sin(relev) * self.dist
    E = np.array((xp, yp, zp))

    self.eye = E
    self.vvec = R - E
    #self.vvec = self.vvec / proj3d.mod(self.vvec)

    if abs(relev) > np.pi/2:
    # upside down                                                                                                                                                                                                                          
      V = np.array((0, 0, -1))
    else:
      V = np.array((0, 0, 1))
    zfront, zback = -self.dist, self.dist

    viewM = proj3d.view_transformation(E, R, V)
    perspM = proj3d.persp_transformation(zfront, zback)
    M0 = np.dot(viewM, worldM)
    M = np.dot(perspM, M0)

    return np.dot(M, scale);

Axes3D.get_proj=get_proj_scale
