#!/usr/bin/env python3


import os
import errno
import copy
import shutil
from time import gmtime, strftime
import numpy as np

def mkdir(thedir):
    try:
        os.makedirs(thedir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

class meshpy2foam():
    ''' Writes large meshes quickly from python to blockMeshDict
    
    Running blockMesh might take a long time then
    A mesh with 860850 vertices took 3.5h to run blockMesh (used 2.8GB of RAM)
    A mesh with 5.9e6 vertices did not finish in 19h (used 17GB of RAM so far)
    A mesh with 33859120 points ran because refinement in x is 334. i.e. it 
    was really 33859120/335=101072 points. blockMesh took close to 2h and 30.5GB, checkMesh took another 6h
    
    blockMesh may say in the beginning:
    Found 1624 undefined faces in mesh; adding to default patch.
    This does not necessarily mean that they will end up in the defaultPatch
    
    bounding box might change a little bit when mirroring and adding the mesh.
    '''
    
    def __init__(self,caseDir):
        # caseDir will be created, the openFoam dir that contanis the mesh
        self.caseDir=caseDir # openFoam case directory
        
        self.v=[] # tuples of vertices will be appended
        self.vCount=0 # none so far
        self.blocks=[]
        self.ref=[]
        self.faces={}
        self.faces['x0']=[] # feel free to rename the boundaries here
        self.faces['x1']=[]
        self.faces['y0']=[]
        self.faces['y1']=[]
        self.faces['z0']=[]
        self.faces['z1']=[]
        self.cells=[]
        self.Prec=6 # precision for coordinates
        
#        print(strftime("%Y-%m-%d_%H:%M:%S", gmtime()) + ' created an instance of meshPython2Foam.')
            
        return
    
    def __add__(self, otherM):
        # overloading + to add two meshes
        bothMeshes = meshpy2foam(self.caseDir)
        nc=len(self.cells)
        for i in range(nc):
            bothMeshes.addCell(self.cells[i][0],self.cells[i][1],self.cells[i][2])
        nc=len(otherM.cells)
        for i in range(nc):
            bothMeshes.addCell(otherM.cells[i][0],otherM.cells[i][1],otherM.cells[i][2])
        return bothMeshes
    
    def addCell(self,v,b,r=(1,1,1)):
        ''' Storing all cell data here before finalize is called'''
        
#        self.cells.append((v,b,r))
        self.cells.append((np.copy(np.round(v,self.Prec)),np.copy(b),np.copy(r)))
        return
    
    def add2Dcell(self,y,z,b,r=(1,1,1)):
        ''' Short cut to add cells that are uniform in x
        
        Before using, define
        self.xmax=
        self.xmin=
        Then supply y and z values of the front view, counter clockwise
        y=np.array([y0,y0,y1,y1])
        z=np.array([z0,z1,z1,z0])
        '''
        
        self.addCell(np.array([
            [[self.xmin],[y[0]],[z[0]]],
            [[self.xmax],[y[0]],[z[0]]],
            [[self.xmax],[y[3]],[z[3]]],
            [[self.xmin],[y[3]],[z[3]]],
            [[self.xmin],[y[1]],[z[1]]],
            [[self.xmax],[y[1]],[z[1]]],
            [[self.xmax],[y[2]],[z[2]]],
            [[self.xmin],[y[2]],[z[2]]]
            ]),b,r)
        return

    def createCell(self,v,b,r=(1,1,1)):
        '''
        Called by finalize, not by the user.
        
        v ... 8 vertices with 3 coordinates each
        b ... boundary strings bounds=['x0','x1','y0','y1','z0','z1'] leave blank '' if not needed
        r ... refinement. Tuple of three integers defining split of that block in every direction
        '''
        
        dx=np.max(v[:,0])-np.min(v[:,0])
        dy=np.max(v[:,1])-np.min(v[:,1])
        dz=np.max(v[:,2])-np.min(v[:,2])
        cellVol=dx*dy*dz
        
        vno=np.zeros(8,dtype=int) # global vertex indeces
        
        if cellVol>0.: # otherwise just completely ignore this cell
            
            # add vertices (every time, will have multiples)
            for i in range(8):
                self.v.append(tuple((v[i,0],v[i,1],v[i,2])))
                vno[i]=self.vCount
                self.vCount +=1 # add in preparation for next
                
            # now add a block
            # numbers in this block refer to vertecies prior to sorting (in finalize)
            self.blocks.append(tuple(vno)) # in the order as passed to addCell()
            self.ref.append(tuple(r)) # refinement
            
            # add faces if applicable:
            if len(b[0])>0: 
                self.faces[b[0]].append(tuple((vno[0],vno[3],vno[7],vno[4])))
            if len(b[1])>0: 
                self.faces[b[1]].append(tuple((vno[1],vno[2],vno[6],vno[5])))
            if len(b[2])>0: 
                self.faces[b[2]].append(tuple((vno[0],vno[1],vno[5],vno[4])))
            if len(b[3])>0: 
                self.faces[b[3]].append(tuple((vno[3],vno[2],vno[6],vno[7])))
            if len(b[4])>0: 
                self.faces[b[4]].append(tuple((vno[0],vno[1],vno[2],vno[3])))
            if len(b[5])>0: 
                self.faces[b[5]].append(tuple((vno[4],vno[5],vno[6],vno[7])))
        
        return 0
    
        
#    def timesNodesBy(self,t):
#        # t=np.array([1.,-1.,1.]) would mirron on y axis
#        # before cells are created, after they are added
#        n=len(self.cells)
#        for i in range(n):
#            self.cells[i][0][:,0] *= t[0]
#            self.cells[i][0][:,1] *= t[1]
#            self.cells[i][0][:,2] *= t[2]
#        return
        
    def addToNodes(self,a):
        # a=np.array([0.,1.,0.]) would shift by 1 along y
        # before cells are created, after they are added
        n=len(self.cells)
        for i in range(n):
            self.cells[i][0][:,0] += np.round(a[0],self.Prec)
            self.cells[i][0][:,1] += np.round(a[1],self.Prec)
            self.cells[i][0][:,2] += np.round(a[2],self.Prec)
        return
    
    def mirrorOnY(self):
        # before cells are created, after they are added
        n=len(self.cells)
        for i in range(n):
            # flip y values around 0.
            self.cells[i][0][:,1] *= -1.
            
            # change order of nodes
            v=np.copy(self.cells[i][0]) # tmp
            self.cells[i][0][0,:]=v[3,:]
            self.cells[i][0][1,:]=v[2,:]
            self.cells[i][0][2,:]=v[1,:]
            self.cells[i][0][3,:]=v[0,:]
            self.cells[i][0][4,:]=v[7,:]
            self.cells[i][0][5,:]=v[6,:]
            self.cells[i][0][6,:]=v[5,:]
            self.cells[i][0][7,:]=v[4,:]
            
            # change boundaries
            y1=self.cells[i][1][2] # may be 'z0' or so
            y0=self.cells[i][1][3]
            if y0=='y1': y0='y0'
            if y1=='y0': y1='y1'
            self.cells[i][1][2]=y0
            self.cells[i][1][3]=y1
        return
    
    def deleteBoundary(self,bs):
        # pass a boundary string, e.g. 'x0'
        n=len(self.cells)
        for i in range(n):
            for j in range(6):
                if self.cells[i][1][j]==bs: self.cells[i][1][j]=''
        return
    
    
    def finalize(self):
        ''' Call this after all cells have been added '''
        
        print(strftime("%Y-%m-%d_%H:%M:%S", gmtime()) + ' beginning meshPython2Foam.finalize() ...')
        
        # mkdir system if not present:
        mkdir(self.caseDir+'/system')
                
        # write aux files necessary to run blockMesh
        self.writeAuxFiles()
        
        # find and print bounding box as info
        xmax=-9e9
        xmin=9e9
        ymax=-9e9
        ymin=9e9
        zmax=-9e9
        zmin=9e9
        nc=len(self.cells)
        for i in range(nc):
            xmin=np.min([xmin,np.min(self.cells[i][0][:,0])])
            xmax=np.max([xmax,np.max(self.cells[i][0][:,0])])
            ymin=np.min([ymin,np.min(self.cells[i][0][:,1])])
            ymax=np.max([ymax,np.max(self.cells[i][0][:,1])])
            zmin=np.min([zmin,np.min(self.cells[i][0][:,2])])
            zmax=np.max([zmax,np.max(self.cells[i][0][:,2])])
        print('Bounding box: [%.8f, %.8f],[%.8f, %.8f],[%.8f, %.8f]' %(xmin,xmax,ymin,ymax,zmin,zmax))
        
        # create cells
        nc=len(self.cells)
        for i in range(nc):
            self.createCell(self.cells[i][0],self.cells[i][1],self.cells[i][2])
        del self.cells # free RAM
        
        # find unique vertices
        self.v,I=np.unique(self.v,return_inverse=True,axis=0)
        
        print(strftime("%Y-%m-%d_%H:%M:%S", gmtime()) + ' sorted mesh.')
        
        # write vertices:
        printstr='(%.'+str(self.Prec)+'f\t%.'+str(self.Prec)+'f\t%.'+str(self.Prec)+'f)\n'
        with open(self.caseDir+'/system/vertices.dat','w') as fid:
            nv=len(self.v)
            for i in range(nv):
                fid.write(printstr %(self.v[i][0],self.v[i][1],self.v[i][2]) )
        print(strftime("%Y-%m-%d_%H:%M:%S", gmtime()) + ' wrote %i vertices.' %nv)
#        del self.v # free RAM
                
        # write blocks:
        with open(self.caseDir+'/system/blocks.dat','w') as fid:
            nb=len(self.blocks)
            for i in range(nb):
                fid.write('hex ( ')
                for j in range(8):
                    fid.write('%i ' %I[self.blocks[i][j]]) # using index after sorting to unique vertices
                fid.write(') ')
                fid.write('( %i %i %i ) ' %( self.ref[i][0],self.ref[i][1],self.ref[i][2] ) )
                fid.write('simpleGrading (1 1 1)\n')
        print(strftime("%Y-%m-%d_%H:%M:%S", gmtime()) + ' wrote %i blocks.' %nb)
#        del self.blocks # free RAM
                
        # write boundary faces:
        for boundary in self.faces:
            nf=len(self.faces[boundary])
            with open(self.caseDir+'/system/faces_'+boundary+'.dat','w') as fid:
                for i in range(nf):
                    fid.write('( ')
                    for j in range(4): # a face has 4 corners (if triangular, two are the same vertex)
                        fid.write('%i ' %(I[self.faces[boundary][i][j]]) )
                    fid.write(')\n')
            print(strftime("%Y-%m-%d_%H:%M:%S", gmtime()) + ' wrote %i faces for boundary ' %nf +boundary+'.')
#        del self.faces # free RAM
        
        
        # for meshVis, create a list of lines at x==xmin
        y0=np.array([])
        y1=np.array([])
        z0=np.array([])
        z1=np.array([])
        nb=len(self.blocks)
        for i in range(nb):
            xbr=self.v[I[self.blocks[i][0]]][0]
            if xbr==xmin:
                # front face is made of up nodes 0,4,7,3. add those 4 lines
                ybr=self.v[I[self.blocks[i][0]]][1] # y bottom right. node 0 in that block
                ytr=self.v[I[self.blocks[i][4]]][1] # y top right. node 4 in that block
                ytl=self.v[I[self.blocks[i][7]]][1] # y top left. node 7 in that block
                ybl=self.v[I[self.blocks[i][3]]][1] # y bottom left. node 3 in that block
                zbr=self.v[I[self.blocks[i][0]]][2] 
                ztr=self.v[I[self.blocks[i][4]]][2] 
                ztl=self.v[I[self.blocks[i][7]]][2] 
                zbl=self.v[I[self.blocks[i][3]]][2]
                
                y0=np.append(y0,ybr); y1=np.append(y1,ytr)
                z0=np.append(z0,zbr); z1=np.append(z1,ztr)
                
                y0=np.append(y0,ytr); y1=np.append(y1,ytl)
                z0=np.append(z0,ztr); z1=np.append(z1,ztl)
                
                y0=np.append(y0,ytl); y1=np.append(y1,ybl)
                z0=np.append(z0,ztl); z1=np.append(z1,zbl)
                
                y0=np.append(y0,ybl); y1=np.append(y1,ybr)
                z0=np.append(z0,zbl); z1=np.append(z1,zbr)
                
        # now delete doublicate lines
        nl=len(y0)
        tmp=np.zeros((nl,4))
        tmp[:,0]=y0
        tmp[:,1]=y1
        tmp[:,2]=z0
        tmp[:,3]=z1
        tmp=np.unique(tmp,axis=0)
        y0=tmp[:,0]
        y1=tmp[:,1]
        z0=tmp[:,2]
        z1=tmp[:,3]
        
        with open(self.caseDir+'/meshVisLines.dat','w') as fid:
            nl=len(y0)
            for i in range(nl):
                fid.write('%.6f\t%.6f\t%.6f\t%.6f\n' %(y0[i], y1[i], z0[i], z1[i]))
        
        
        print(strftime("%Y-%m-%d_%H:%M:%S", gmtime()) + ' finished meshPython2Foam.finalize() ...')
        
        return 0
    
    def writeAuxFiles(self):
        
        thisDir=os.path.dirname(__file__)
        shutil.copy(thisDir+'/blockMeshDict', self.caseDir+'/system/')
        shutil.copy(thisDir+'/controlDict',   self.caseDir+'/system/')
        shutil.copy(thisDir+'/fvSchemes',     self.caseDir+'/system/')
        shutil.copy(thisDir+'/fvSolution',    self.caseDir+'/system/')
        
        return 0
    
def half2N_y0(m,ndoub):
    '''
    Half a mesh to n meshes: add a copy of m mirrored on y0, 
    then doublicate the combined mesh ndoub times in y.
    '''
    
    # find max(y)
    ymax=0.
    n=len(m.cells)
    for i in range(n):
        ymax=np.max([ymax,np.max(m.cells[i][0][:,1])])
    ymax=np.round(ymax,6)
    
    # create a copy and flip it
    m2=copy.deepcopy(m) # otherwise changes to m2 affect m
    m2.mirrorOnY()
    m2.deleteBoundary('y1')
    m1=copy.deepcopy(m)
    m1.deleteBoundary('y0')
    m=m1+m2
    m.addToNodes(np.array([0.,ymax,0.]))
    
    # multiply in y:
    m1=copy.deepcopy(m)
    m1.deleteBoundary('y0')
    for r in range(ndoub-1):
        m11=copy.deepcopy(m1)
        m11.addToNodes(np.array([0.,(r+1)*2.*ymax,0.]))
        m.deleteBoundary('y1')
        m += m11
    
    return m

def half2N_y1(m,ndoub):
    '''
    Half a mesh to n meshes: add a copy of m mirrored on y1, 
    then doublicate the combined mesh ndoub times in y.
    '''
    
    # find max(y)
    ymax=0.
    n=len(m.cells)
    for i in range(n):
        ymax=np.max([ymax,np.max(m.cells[i][0][:,1])])
    ymax=np.round(ymax,6)
    
    # create a copy and flip it
    m2=copy.deepcopy(m) # otherwise changes to m2 affect m
    m2.mirrorOnY()
    m2.deleteBoundary('y0')
    m1=copy.deepcopy(m)
    m1.deleteBoundary('y1')
    m2.addToNodes(np.array([0.,2*ymax,0.]))
    m=m1+m2
    
    # multiply in y:
    m1=copy.deepcopy(m)
    m1.deleteBoundary('y0')
    for r in range(ndoub-1):
        m11=copy.deepcopy(m1)
        m11.addToNodes(np.array([0.,(r+1)*2.*ymax,0.]))
        m.deleteBoundary('y1')
        m += m11
        
    return m