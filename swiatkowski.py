import numpy as np
import sympy as sym
import itertools as it
import collections as col
import copy
from operator import itemgetter
import os
import time
#from sage.all import *

#################INPUT############
#a graph without self-loops
#number of particles
Npart=2
#n-th homology
nhom=1

#################################################################
#################################################################
'''
#D4
name='D4'
N0=4
graph=np.zeros((N0,N0),dtype=int)
graph[0][1]=2
graph[0][3]=2
graph[1][2]=2
graph[2][3]=2
'''
'''
#crate4x4
name='crate4x4'
N0=12
graph=np.zeros((N0,N0),dtype=int)
graph[0][1]=1
graph[0][7]=1
graph[0][8]=1
graph[1][2]=1
graph[1][8]=1
graph[2][3]=1
graph[2][9]=1
graph[3][4]=1
graph[3][9]=1
graph[4][5]=1
graph[4][10]=1
graph[5][6]=1
graph[5][10]=1
graph[6][11]=1
graph[6][7]=1
graph[7][11]=1
graph[8][9]=1
graph[8][11]=1
graph[9][10]=1
graph[10][11]=1

'''
'''
#K2,p
name='K2,4'
N0=2
graph=np.zeros((N0,N0),dtype=int)
graph[0][1]=7
'''
'''
#lasso
name='lasso'
N0=3
graph=np.zeros((N0,N0),dtype=int)
graph[0][1]=1
graph[1][2]=2
'''
'''
#K4
name='K4'
print 'K4'
N0=4
graph=np.zeros((N0,N0),dtype=int)
graph[0][1]=1
graph[0][2]=1
graph[0][3]=1
graph[1][2]=1
graph[1][3]=1
graph[2][3]=1
'''
'''
#K5
print 'K5'
name="K5"
N0=5
graph=np.zeros((N0,N0),dtype=np.int8)
graph[0][1]=1
graph[0][2]=1
#graph[0][3]=1
#graph[0][4]=1
graph[1][2]=1
graph[1][3]=1
graph[1][4]=1
graph[2][3]=1
graph[2][4]=1
graph[3][4]=1
'''
'''
#K5-K5
name='K5-K5'
print name
N0=10
graph=np.zeros((N0,N0),dtype=np.int8)
graph[0][1]=1
graph[0][2]=1
graph[0][3]=1
graph[0][4]=1
graph[1][2]=1
graph[1][3]=1
graph[1][4]=1
graph[2][3]=1
graph[2][4]=1
graph[3][4]=1
graph[0][5]=1
graph[5][6]=1
graph[5][7]=1
graph[5][8]=1
graph[5][9]=1
graph[6][7]=1
graph[6][8]=1
graph[6][9]=1
graph[7][8]=1
graph[7][9]=1
graph[8][9]=1
'''
'''
#K3,3
N0=6
name='K3,3'
graph=np.zeros((N0,N0),dtype=np.int8)
graph[0][1]=1
graph[0][3]=1
graph[0][4]=1
graph[1][2]=1
graph[1][5]=1
graph[2][3]=1
graph[2][4]=1
graph[3][5]=1
graph[4][5]=1
'''
'''
#K6
print 'K6'
name="K6"
N0=6
graph=np.zeros((N0,N0),dtype=np.int8)
graph[0][1]=1
graph[0][2]=1
graph[0][3]=1
graph[0][4]=1
graph[0][5]=1
graph[1][2]=1
graph[1][3]=1
graph[1][4]=1
graph[1][5]=1
graph[2][3]=1
graph[2][4]=1
graph[2][5]=1
graph[3][4]=1
graph[3][5]=1
graph[4][5]=1
'''
'''
#petersen1
print 'petersen1'
name='petersen1'
N0=7
graph=np.zeros((N0,N0),dtype=np.int8)
graph[0][1]=1
graph[0][3]=1
graph[0][4]=1
graph[0][5]=1
graph[0][6]=1
graph[1][2]=1
graph[1][4]=1
graph[1][6]=1
graph[2][3]=1
graph[2][5]=1
graph[3][4]=1
graph[3][6]=1
graph[4][5]=1
graph[4][6]=1
graph[5][6]=1
'''
#Y
print 'Y'
name='Y'
N0=4
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[1][2]=1
graph[1][3]=1

################################
graph=graph+np.transpose(graph)
#################################################################
#################################################################

def edge_enum():
    edges=[[[] for i in xrange(N0)] for j in range(N0)]
    ne=0
    for i in xrange(N0):
        for j in xrange(i,N0,1):
            if graph[i][j]>0:
                edges[i][j]+=range(ne,ne+graph[i][j],1)
                ne+=graph[i][j]
        for j in range(0,i,1):
            if graph[i][j]>0:
                edges[i][j]+=copy.deepcopy(edges[j][i])
                #edges[i][j]+=copy(edges[j][i])
    return (edges,ne)
    

def generators(d):
    hedges=[]
    (edges_matrix,ne)=edge_enum()
    dcells=[]
    for v in xrange(N0):
        edges_v=np.sum(edges_matrix[v])
        #hedges.append([(evpair,v) for evpair in it.combinations(np.sum(edges_matrix[v]),2)])
        hedges.append([((edges_v[0],e),v) for e in edges_v[1:]])
    dcells=[list(c) for vsubset in it.combinations(xrange(N0),d) for c in it.product(*[hedges[nv] for nv in vsubset])]
    free_part_distr=[col.Counter(com) for com in it.combinations_with_replacement(xrange(ne),Npart-d)]
    return it.product(dcells,free_part_distr)
    
def boundary(cell,d):
    sign=1
    boundary=[]
    for ih in xrange(d):
        distin=copy.deepcopy(cell[1])
        #distin=copy(cell[1])
        distin.update([cell[0][ih][0][0]])
        distout=copy.deepcopy(cell[1])
        #distout=copy(cell[1])
        distout.update([cell[0][ih][0][1]])
        boundary.append(((cell[0][:ih]+cell[0][ih+1:],distin),sign))
        boundary.append(((cell[0][:ih]+cell[0][ih+1:],distout),-sign))
        sign*=(-1)
    return boundary
    
def boundary_operator(d):
    genminus=list(generators(d-1))
    #print len(genminus)
    mtx=[]
    for c in generators(d):
        bnd=boundary(c,d)
        mtx.append(np.zeros(len(genminus),dtype=np.int8))
        for bndc in bnd:
            mtx[-1][genminus.index(bndc[0])]+=bndc[1] 
    return np.array(mtx)

def find_cycles(dimension):
    cellsd=list(generators(dimension))
    bnd_op=boundary_operator(dimension)
    cycles=[]
    bndMat=sym.Matrix(np.transpose(bnd_op))
    kernel=bndMat.nullspace()
    for cycle in kernel[0:]:
        cells=[]
        for i in xrange(len(cycle)):
            if cycle[i]!=0:
                cells.append((cellsd[i],int(cycle[i])))
        cycles.append(cells)
    return cycles    

start=time.time()

kermtx=boundary_operator(nhom)
print np.shape(kermtx)
print kermtx.nbytes
dimkerDn=np.shape(kermtx)[0]-np.linalg.matrix_rank(kermtx)
print "dimkerDn="+str(dimkerDn)
'''
for cycle in find_cycles(nhom):
    for cell in cycle: print cell
'''
if Npart-(nhom+1)>=0:
    immtx=boundary_operator(nhom+1)
    print np.shape(immtx)
    print immtx.nbytes
    dimImDnplus1=np.linalg.matrix_rank(immtx)
    print "dimImDnplus1="+str(dimImDnplus1)
    print 'betti number='+str(dimkerDn-dimImDnplus1)
    '''
    #smith normal form in sage
    Aim=matrix(ZZ,immtx,sparse=True)
    Aimsmith=(Aim.smith_form())[0]
    smithdiag=np.array((Aimsmith).diagonal())
    '''
#print 'torsion='+str([a for a in smithdiag if a >1])

print 'elapsed time: '+str(int(time.time()-start))+'s'
    
np.savetxt(os.path.expanduser('~/Dropbox/quantum graphs/numerics/Dimage'+str(name)+'_h'+str(nhom)+'_n'+str(Npart)+'.txt'),immtx,fmt='%d',delimiter='\t')


