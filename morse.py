import numpy as np
import sympy as sym
import itertools as it
import collections as col
import copy
from operator import itemgetter
import os


#################
#number of particles
Npart=7
#n-th homology
nhom=5

#################################################################
#################################################################
'''
#bowtie
N0=5
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][2]=-1
graph[1][2]=1
graph[1][3]=1
graph[1][4]=1
graph[3][4]=-1
'''

#K3,3
N0=6
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][3]=-1
graph[0][4]=-1
graph[1][2]=1
graph[1][5]=-1
graph[2][3]=1
graph[2][4]=1
graph[3][5]=-1
graph[4][5]=1

'''
#K5
print 'K5'
N0=5
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][2]=-1
graph[0][3]=-1
graph[0][4]=-1
graph[1][2]=1
graph[1][3]=-1
graph[1][4]=-1
graph[2][3]=1
graph[2][4]=-1
graph[3][4]=1
'''
'''
#K4
print 'K4'
N0=4
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][2]=-1
graph[0][3]=-1
graph[1][2]=1
graph[1][3]=-1
graph[2][3]=1
'''
'''
#myexample
N0=7
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][2]=-1
graph[0][4]=-1
graph[1][2]=1
graph[1][3]=1
graph[1][5]=-1
graph[1][6]=-1
graph[2][3]=-1
graph[2][4]=-1
graph[3][6]=-1
graph[3][4]=1
graph[4][5]=1
graph[4][6]=1
graph[5][6]=-1
'''
'''
#K3,3_modified
N0=6
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][3]=-1
graph[1][2]=1
graph[2][3]=1
graph[2][4]=1
graph[3][5]=-1
graph[4][5]=1
'''
'''
#myexample_modified
N0=7
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][2]=-1
graph[0][4]=-1
graph[1][2]=1
graph[1][3]=1
graph[1][5]=-1
graph[1][6]=-1
graph[2][4]=-1
graph[3][4]=1
graph[4][5]=1
graph[4][6]=1
graph[5][6]=-1
'''

graph=graph+np.transpose(graph)

#################################################################
#################################################################    
def expand_graph(graph0,e_to_expand,val):
    graph_new=copy.deepcopy(graph0)
    graph_new[e_to_expand[0]][e_to_expand[1]]=0
    graph_new[e_to_expand[1]][e_to_expand[0]]=0
    if graph0[e_to_expand[0]][e_to_expand[1]]>0:
        new_v=e_to_expand[1]+1
        graph_new=np.insert(graph_new,e_to_expand[1],0,axis=1)
        graph_new=np.insert(graph_new,e_to_expand[1],0,axis=0)
        graph_new[e_to_expand[0]][e_to_expand[1]]=val
        graph_new[e_to_expand[1]][e_to_expand[0]]=val
    else:
        outgoing=[]
        is_branch=False
        for iv in np.arange(e_to_expand[1],len(graph0[:,e_to_expand[1]]),1):
            if graph_new[iv][e_to_expand[1]]==1:
                outgoing.append(iv)
                is_branch=True
        if is_branch:
            new_v=min(outgoing)
        else:
            new_v=e_to_expand[1]+1
        graph_new=np.insert(graph_new,new_v,0,axis=1)
        graph_new=np.insert(graph_new,new_v,0,axis=0)
        graph_new[e_to_expand[0]][new_v]=-val
        graph_new[new_v][e_to_expand[0]]=-val
    graph_new[e_to_expand[1]][new_v]=1
    graph_new[new_v][e_to_expand[1]]=1
    return graph_new

def divide_edges(graph):
    size=graph.shape[0]
    tree=[]
    del_edges=[]
    blocked_edges=[]
    vertices=[i+1 for i in xrange(size)]
    for com in it.combinations(xrange(size),2):
        if graph[com[0]][com[1]]>0:
            tree.append((com[0]+1,com[1]+1))
        elif graph[com[0]][com[1]]<0:
            del_edges.append((com[0]+1,com[1]+1))
    for i in xrange(size):
        if (graph[i]>0).sum()>2:
            blocked=[]
            for e in tree:
                if e[0]==vertices[i]:
                    blocked.append(e)
            for comb in it.combinations(blocked,2):
                blocked_edges.append((comb[1],comb[0][1]))
    return (tree,blocked_edges,del_edges,vertices)
    
   
def subdivision(Npart,graph0):
    graph=copy.deepcopy(graph0)
    
    for i in xrange(N0):
        for j in range(i+1,N0,1):
            if graph0[i][j]!=0 and (graph0[j] != 0).sum()!=2 and (graph0[i] != 0).sum()!=2:
                    graph[i][j]*=2
                    graph[j][i]*=2
            elif (graph0[i] != 0).sum()==2 and graph0[i][j]==-1:
                rowpart=list(graph0[i])
                graph[min(rowpart.index(1),i)][max(rowpart.index(1),i)]*=2
                graph[max(rowpart.index(1),i)][min(rowpart.index(1),i)]*=2
    n_expansions=1
    while True:
        size=graph.shape[0]
        expanded=False
        for com in it.combinations(xrange(size),2):
            if abs(graph[com[0]][com[1]])==2:
                    if Npart>2:
                        if n_expansions%(Npart-2)==0:
                            val=1
                        else:
                            val=2
                    else:
                        continue
                    graph=expand_graph(graph,com,val)
                    expanded=True
                    n_expansions+=1
                    break
            if expanded:
                break
        if expanded==False:
            break
    return divide_edges(graph)
    

(tree,blocked_edges,del_edges,vert)=subdivision(Npart,graph)
print 'tree: \t'+str(tree)
print 'deleted edges: \t'+str(del_edges)
########################################################

def find_crit_cells(dim):
    cells=set()
    n_del_min=2*dim-Npart
    if n_del_min<0:
        n_del_min=0
    for n_del in range(n_del_min,dim+1,1):
        del_coll=it.combinations(del_edges,n_del)
        blocked_coll=it.combinations(blocked_edges,dim-n_del)
        for elem in it.product(del_coll,blocked_coll):
            flat=[]
            for i in xrange(n_del):
                flat.append(elem[0][i][0])
                flat.append(elem[0][i][1])
            for i in xrange(dim-n_del):
                flat.append(elem[1][i][0][0])
                flat.append(elem[1][i][0][1])
                flat.append(elem[1][i][1])
            if len(set(flat))==len(flat):
                for comb_verts in it.combinations(tuple(set(vert)-set(flat)),Npart+n_del-2*dim):
                    all_blocked=True
                    for v in comb_verts:
                        if v==1 or (next((e[0] for e in tree if e[1]==v),None) in (flat+list(comb_verts))):
                            continue
                        all_blocked=False
                        break
                    if all_blocked:
                        candidate=elem[0]
                        candidate_v=comb_verts
                        for item in elem[1]:
                            candidate+=(item[0],)
                            candidate_v+=(item[1],)
                        cells.add(sorted_cell(candidate+candidate_v,dim))
    return cells

    
def boundary(cell,dim):
    boundary_cells={}
    sign=-1
    for i in xrange(dim):
        boundary_cells.update({cell[:i]+cell[i+1:dim]+tuple(sorted(cell[dim:]+(cell[i][1],))):sign})
        boundary_cells.update({cell[:i]+cell[i+1:dim]+tuple(sorted(cell[dim:]+(cell[i][0],))):-sign})
        sign*=-1
    return boundary_cells


def sorted_cell(cell, dim):
    return tuple(sorted(cell[:dim], key=itemgetter(0,1)))+tuple(sorted(cell[dim:]))    
    
def no_order_resp_edge(cell,dim):
    edges_ordered=sorted(cell[:dim], key=itemgetter(1,0))
    for ic in xrange(dim):
        if (edges_ordered[ic] in del_edges):
            continue
        isblocked=False
        for cb in blocked_edges:
            if cb[0]==edges_ordered[ic] and cb[1] in cell[dim:]:
                isblocked=True
                break
        if isblocked==False:
            return (False,edges_ordered[ic]) #returns the minimal o.r. edge
    return (True,'x')

def all_vert_blocked(cell,dim):
    flat=tuple([item for sublist in cell[:dim] for item in sublist])+cell[dim:]
    for v in cell[dim:]:
        if v==1 or next((e[0] for e in tree if e[1]==v),None) in flat:
            continue 
        return (False,v) #returns smallest unblocked vertex
    return (True,'x')
    

def iscritical(cell, dim):
    if no_order_resp_edge(cell,dim)[0] and all_vert_blocked(cell,dim)[0]:
        return True
    else:
        return False
        
def principal_reduction(cell,dim,v):
    #print v
    reduced=cell[:dim]
    reduced+=(next(e for e in tree if e[1]==v),)
    reduced=sorted_cell(reduced,dim+1)
    ind=cell.index(v)
    reduced+=cell[dim:ind]+cell[ind+1:]
    return reduced
            
def boundary_m(cell,dim,criticals):
    bnd_old=copy.deepcopy(boundary(cell,dim))
    while True:
        bnd_new={}
        for c in bnd_old:
            if c in criticals:
                if c in bnd_new:
                    bnd_new[c]+=bnd_old[c]
                else:
                    bnd_new[c]=bnd_old[c]
            elif no_order_resp_edge(c,dim-1)[0]:
                unblocked=all_vert_blocked(c,dim-1)
                if unblocked[0]==False:
                    step=principal_reduction(c,dim-1,unblocked[1])
                    bnd=boundary(step,dim)
                    for cb in bnd:
                        if cb!=c:
                            if cb in bnd_new:
                                bnd_new[cb]+=(-bnd[c])*bnd_old[c]*bnd[cb]
                            else:
                                bnd_new[cb]=(-bnd[c])*bnd_old[c]*bnd[cb]
            else:
                e=no_order_resp_edge(c,dim-1)[1]
                flat=tuple([item for sublist in c[:dim-1] for item in sublist])+c[dim-1:]
                for v in c[dim-1:]:
                    if v!=1 and v<e[1] and (next((e[0] for e in tree if e[1]==v),None) in flat)==False:
                        step=principal_reduction(c,dim-1,v)
                        bnd=boundary(step,dim)
                        for cb in bnd:
                            if cb!=c:
                                if cb in bnd_new:
                                    bnd_new[cb]+=-bnd[c]*bnd_old[c]*bnd[cb]
                                else:
                                    bnd_new[cb]=-bnd[c]*bnd_old[c]*bnd[cb]
                        break
        all_crit=True
        for c in bnd_new:
            if (c in criticals)==False:
                all_crit=False
                break
        if all_crit:
            return bnd_new
        else:
            bnd_old=copy.deepcopy(bnd_new)
            
                                
def boundary_operator(dim): #dim-dimension of the domain
    cells_dom=find_crit_cells(dim)
    cells_im_set=find_crit_cells(dim-1)
    cells_im={}
    ind_im=0
    for c in cells_im_set:
         cells_im.update({c:ind_im})
         ind_im+=1
    Dmtx=np.zeros((len(cells_dom),len(cells_im)))
    ind_d=0
    n_dom=len(cells_dom)
    for cd in cells_dom:
        if (int(100.*float(ind_d+1)/n_dom)-int(100.*float(ind_d)/n_dom))==1:
            print 'boundary_operator:\t'+str(int(100.*float(ind_d)/n_dom))+'proc.'
        bnd=boundary_m(cd,dim,cells_im_set)
        for cb in bnd.items():
            Dmtx[ind_d][cells_im[cb[0]]]=cb[1]
        ind_d+=1
    return Dmtx
    
def find_cycles(bnd_op,cells_dom):
    cellsd=[]
    cycles=[]
    for cd in cells_dom:
        cellsd.append(cd)
    bndMat=sym.Matrix(np.transpose(bnd_op))
    kernel=bndMat.nullspace()
    for cycle in kernel[0:]:
        cells={}
        for i in xrange(len(cycle)):
            if cycle[i]!=0:
                cells.update({cellsd[i]:int(cycle[i])})
        cycles.append(cells)
    return cycles
        
        
######################################################
cellsn=find_crit_cells(nhom)
print "no. of critical "+str(nhom)+"-cells "+str(len(cellsn))
cellsnminus=find_crit_cells(nhom-1)
print "no. of critical "+str(nhom-1)+"-cells "+str(len(cellsnminus))
cellsnplus=find_crit_cells(nhom+1)
print "no. of critical "+str(nhom+1)+"-cells "+str(len(cellsnplus))

Dnmtx=boundary_operator(nhom)
print "rank of boundary "+str(np.linalg.matrix_rank(Dnmtx))                    
print "dimkerDn="+str(len(cellsn)-np.linalg.matrix_rank(Dnmtx))

n_cycles=find_cycles(Dnmtx,cellsn)
f_cycles=open(os.path.expanduser('~/Dropbox/quantum graphs/numerics/cycles.txt'),'w')
for c in n_cycles:
    f_cycles.write(str(c)+'\n\n')
    print c
f_cycles.close()

f_nplus=open(os.path.expanduser('~/Dropbox/quantum graphs/numerics/nplus_cells.txt'),'w')
for c in cellsnplus:
    f_nplus.write(str(c)+'\n')
    f_nplus.write(str(boundary_m(c,nhom+1,cellsn))+'\n\n')
f_nplus.close()

Dnpmtx=boundary_operator(nhom+1)              

to_delete=[]
flag=True
n_deleted=0

while flag:
    for i in xrange(Dnpmtx.shape[1]):
        if len(col.Counter(Dnpmtx[:,i]))==1 and Dnpmtx[0][i]==0:
            Dnpmtx=np.delete(Dnpmtx,i,1)
            n_deleted+=1
            break
        flag=False

print "size of Dimage "+str(Dnpmtx.shape)
print "deleted columns: "+str(n_deleted)

np.savetxt(os.path.expanduser('~/Dropbox/CFT/quantum graphs/numerics/Dimage.txt'),Dnpmtx,fmt='%d',delimiter='\t')

