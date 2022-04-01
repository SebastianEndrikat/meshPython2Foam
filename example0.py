#!/usr/bin/env python3

'''
Create a mesh of some silly arch to demonstrate how to add 
cells to the grid and assign boundary faces
'''

import sys
sys.path.insert(0, './src')
import meshPython2Foam
import numpy as np


# create instance of the mesh building class:
foamDir='example0'
m=meshPython2Foam.meshpy2foam(foamDir)


# mesh settings
ref=(10,1,1) # refinement in blockMesh
x0=0.
x1=2. # extend in x happens to be the same for all cells in the mesh


################# add a bunch of cells in one area
y0=0.; y1=1.
bounds=['x0','x1','y0','z0','',''] # boundary zone names for every block. x0,x1,y0,y1,z0,z1
# bounds[3] is the bottom boundary for this block even tho it's on the left side of the block
z=np.array([0,1.5,2,3,5]) # may be non-uniform
for k in range(0,len(z)-1):
    z0=z[k]
    z1=z[k+1]
    if k==0: 
        bounds[4]='z0'
    else:
        bounds[4]='' # inside the mesh, not a boundary
    m.addCell(np.array([ 
                 [x0,y0,z0],
                 [x1,y0,z0],
                 [x1,y1,z0],
                 [x0,y1,z0],
                 [x0,y0,z1],
                 [x1,y0,z1],
                 [x1,y1,z1],
                 [x0,y1,z1]]),bounds,ref)

    
    
################# add a bunch of cells in another area
y0=3.; y1=4.
bounds=['x0','x1','z0','y1','',''] # boundary zone names for every block. x0,x1,y0,y1,z0,z1
# bounds[2] is the bottom boundary for this block even tho it's on the right side of the block
for k in range(0,len(z)-1):
    z0=z[k]
    z1=z[k+1]
    if k==0: 
        bounds[4]='z0'
    else:
        bounds[4]='' # inside the mesh, not a boundary
    m.addCell(np.array([ 
                 [x0,y0,z0],
                 [x1,y0,z0],
                 [x1,y1,z0],
                 [x0,y1,z0],
                 [x0,y0,z1],
                 [x1,y0,z1],
                 [x1,y1,z1],
                 [x0,y1,z1]]),bounds,ref)
    
    
################# add a bunch of cells on top of the previous two sections
y=np.linspace(0,4,5)
z0=np.max(z)
z1=z1+1.
for j in range(0,len(y)-1):
    y0=y[j]
    y1=y[j+1]
    bounds=['x0','x1','','','','z1'] # boundary zone names for every block. x0,x1,y0,y1,z0,z1
    if j==0: 
        bounds[2]='y0'
    if y1==np.max(y):
        bounds[3]='y1'
    if y0>=1 and y1<=3:
        bounds[4]='z0'
    m.addCell(np.array([ 
                 [x0,y0,z0],
                 [x1,y0,z0],
                 [x1,y1,z0],
                 [x0,y1,z0],
                 [x0,y0,z1],
                 [x1,y0,z1],
                 [x1,y1,z1],
                 [x0,y1,z1]]),bounds,ref)
    

m.finalize() # write foam files to disk
print('Done.')
print('Now go and cd '+foamDir+'; blockMesh; checkMesh; paraFoam; '+
      '(and optionally also foamMeshToFluent to convert the whole mesh to a fluent mesh)')
