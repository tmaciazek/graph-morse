import numpy as np
import sympy as sym
import itertools as it
import collections as col
import copy
from operator import itemgetter
import os
import time


#################
#number of particles
Npart=2
#n-th homology
nhom=1

#################################################################
#################################################################
'''
#K3,3
name='K3,3'
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
'''
#K5
print 'K5'
name="K5"
N0=6
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
#K5
print 'K5'
name="K5"
N0=8
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][7]=-1
graph[1][2]=1
graph[1][5]=-1
graph[1][6]=-1
graph[2][3]=1
graph[2][4]=1
graph[3][7]=-1
graph[4][6]=-1
graph[2][5]=1
graph[5][7]=-1
graph[5][6]=1
graph[6][7]=1
'''
#K4
print 'K4'
N0=4
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][2]=-1
graph[0][3]=-1
graph[1][2]=1
graph[1][3]=1
graph[2][3]=-1
'''

'''
#K4
print 'K4'
N0=7
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][5]=-1
graph[0][6]=-1
graph[1][2]=1
graph[1][4]=-1
graph[2][3]=1
graph[2][6]=1
graph[3][4]=1
graph[3][5]=1
'''
'''
#wheel4
name='wheel4'
print name
N0=5
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][2]=-1
graph[0][4]=-1
graph[1][2]=1
graph[1][3]=1
graph[1][4]=1
graph[2][3]=-1
graph[3][4]=-1
'''
'''
#wheel5
print 'wheel5'
name='wheel5'
N0=6
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][4]=-1
graph[0][5]=-1
graph[1][2]=1
graph[1][5]=-1
graph[2][3]=1
graph[2][5]=-1
graph[3][4]=1
graph[3][5]=-1
graph[4][5]=1
'''
'''
#F2,3 aka lasso
N0=4
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[1][2]=1
graph[1][3]=1
graph[2][3]=-1
'''
'''
#baloon
N0=8
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[1][2]=1
graph[1][5]=1
graph[2][3]=1
graph[2][4]=1
graph[3][4]=-1
graph[5][6]=1
graph[5][7]=1
graph[6][7]=-1
'''
'''
#triple_torus or K2,4
N0=6
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][3]=-1
graph[0][4]=-1
graph[0][5]=-1
graph[1][2]=1
graph[2][3]=1
graph[2][4]=1
graph[2][5]=1
'''
'''
#theta-4
print 'Theta'
N0=8
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][5]=-1
graph[0][7]=-1
graph[1][2]=1
graph[2][3]=1
graph[3][4]=1
graph[3][6]=1
graph[4][5]=1
graph[6][7]=1
'''
'''
#theta-5
print 'Theta'
N0=11
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][7]=-1
graph[0][10]=-1
graph[1][2]=1
graph[2][3]=1
graph[3][4]=1
graph[4][5]=1
graph[4][8]=1
graph[5][6]=1
graph[6][7]=1
graph[8][9]=1
graph[9][10]=1
'''
'''
#theta4-3
print 'Theta'
N0=10
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][5]=-1
graph[0][7]=-1
graph[0][9]=-1
graph[1][2]=1
graph[2][3]=1
graph[3][4]=1
graph[3][6]=1
graph[3][8]=1
graph[4][5]=1
graph[6][7]=1
graph[8][9]=1
'''
##############################
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

    
def boundary2(cell):
    letters=[None]*4
    exps=[None]*4
    letters[3]=(cell[:0]+cell[1:2]+tuple(sorted(cell[2:]+(cell[0][1],))))
    exps[3]=1
    letters[1]=(cell[:0]+cell[1:2]+tuple(sorted(cell[2:]+(cell[0][0],))))
    exps[1]=-1
    letters[2]=(cell[:1]+cell[2:2]+tuple(sorted(cell[2:]+(cell[1][1],))))
    exps[2]=-1
    letters[0]=(cell[:1]+cell[2:2]+tuple(sorted(cell[2:]+(cell[1][0],))))
    exps[0]=1
    return (letters,exps)


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
    return (True, None)
    

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
            
def reduction2(letters,exps,criticals1):
    letters_new=[]
    exps_new=[]
    for i in xrange(len(letters)):
        if letters[i] in criticals1:
            letters_new.append(letters[i])
            exps_new.append(exps[i])
        elif no_order_resp_edge(letters[i],1)[0]:
            unblocked=all_vert_blocked(letters[i],1)
            if unblocked[0]==False:
                match=principal_reduction(letters[i],1,unblocked[1])
                bnd=boundary2(match)
                ib=bnd[0].index(letters[i])
                if bnd[1][ib]>0:
                    lett_toinsert=[bnd[0][(ib-j-1)%4] for j in xrange(3)]
                    exps_toinsert=[-bnd[1][(ib-j-1)%4] for j in xrange(3)]
                else:
                    lett_toinsert=[bnd[0][(ib+j+1)%4] for j in xrange(3)]
                    exps_toinsert=[bnd[1][(ib+j+1)%4] for j in xrange(3)]
                if exps[i]>0:
                    letters_new+=lett_toinsert
                    exps_new+=exps_toinsert
                else:
                    lett_toinsert.reverse()
                    exps_toinsert.reverse()
                    exps_toinsert=[-sgn for sgn in exps_toinsert]
                    letters_new+=lett_toinsert
                    exps_new+=exps_toinsert        
        else:
            e=no_order_resp_edge(letters[i],1)[1]
            flat=tuple([item for sublist in (letters[i])[:1] for item in sublist])+(letters[i])[1:]
            for v in letters[i][1:]:
                if v!=1 and v<e[1] and (next((e[0] for e in tree if e[1]==v),None) in flat)==False:
                    match=principal_reduction(letters[i],1,v)
                    bnd=boundary2(match)
                    ib=bnd[0].index(letters[i])
                    if bnd[1][ib]>0:
                        lett_toinsert=[bnd[0][(ib-j-1)%4] for j in xrange(3)]
                        exps_toinsert=[-bnd[1][(ib-j-1)%4] for j in xrange(3)]
                    else:
                        lett_toinsert=[bnd[0][(ib+j+1)%4] for j in xrange(3)]
                        exps_toinsert=[bnd[1][(ib+j+1)%4] for j in xrange(3)]
                    if exps[i]>0:
                        letters_new+=lett_toinsert
                        exps_new+=exps_toinsert
                    else:
                        lett_toinsert.reverse()
                        exps_toinsert.reverse()
                        letters_new+=lett_toinsert
                        exps_new+=[-sgn for sgn in exps_toinsert]
                    break
    return (letters_new,exps_new)
        

######################################################
cellsnminus=find_crit_cells(0)
print "no. of critical "+str(0)+"-cells "+str(len(cellsnminus))
cells1=find_crit_cells(1)
gens=list(cells1)
print "no. of generators: "+str(len(cells1))
cells2=find_crit_cells(2)
print "no. of relators: "+str(len(cells2))

for c2 in cells2:
    print c2
    bndc=boundary2(c2)
    letti,expsi=bndc[0],bndc[1]
    while True:
        red=reduction2(letti,expsi,gens)
        letti,expsi=red[0],red[1]
        if all([el in gens for el in letti]):
            break
    print [('g'+str(gens.index(letti[i])),expsi[i]) for i in xrange(len(letti))]

for i in xrange(len(gens)):
    print 'g'+str(i)+':'+str(gens[i])
'''
#gamma
#letti=[((5,9),1,2,4),((5, 6), 1, 2, 4),((5, 9), 1, 2, 6)]
#letti=[((5,9),4),((5, 6),4),((5, 4),6),((5, 9),6),((5, 6),9),((5, 4),9)]
#letti=[((5,9),4,7,8),((5, 6),4,7,8),((5, 4),6,7,8),((5, 9),6,7,8),((5, 6),7,8,9),((5, 4),7,8,9)]
#letti=[((5,9),1,4,7),((5, 6),1,4,7),((5, 4),1,6,7),((5, 9),1,6,7),((5, 6),1,7,9),((5, 4),1,7,9)]
#letti=[((5,9),1,4,10),((5, 6),1,4,10),((5, 4),1,6,10),((5, 9),1,6,10),((5, 6),1,9,10),((5, 4),1,9,10)]
#letti=[((5,9),4,10,11),((5, 6),4,10,11),((5, 4),6,10,11),((5, 9),6,10,11),((5, 6),9,10,11),((5, 4),9,10,11)]
#letti=[((5,9),4,10),((5, 6),4,10),((4, 5),6,10),((5, 9),6,10),((5, 6),9,10),((4, 5),9,10)]
#expsi=[1,-1,-1,-1,1,1]

#letti=[((1,11),2),((1, 8),2),((1, 2),8),((1, 11),8),((1, 8),11),((1, 2),11)]
#expsi=[1,-1,1,-1,1,-1]


#alphaU
#letti=[((1,8),9,10)]
letti=[((1,2),9,10),((2,3),9,10),((3,4),9,10),((4,5),9,10),((5,6),9,10),((6,7),9,10),((7,8),9,10),((1,8),9,10)]
#letti=[((1,2),9,10,11),((2,3),9,10,11),((3,4),9,10,11),((4,5),9,10,11),((5,6),9,10,11),((6,7),9,10,11),((7,8),9,10,11),((1,8),9,10,11)]
expsi=[1,1,1,1,1,1,1,-1]
'''

'''
#alphaD
#letti=[((5,6),2,3,4),((6,7),2,3,4),((7,8),2,3,4),((1,8),2,3,4),((1,11),2,3,4),((10,11),2,3,4),((9,10),2,3,4),((5,9),2,3,4)]
letti=[((5,6),2,3),((6,7),2,3),((7,8),2,3),((1,8),2,3),((1,11),2,3),((10,11),2,3),((9,10),2,3),((5,9),2,3)]
expsi=[-1,-1,-1,1,-1,1,1,1]
'''
'''
#alphaO
letti=[((1,2),6,7,8),((2,3),6,7,8),((3,4),6,7,8),((4,5),6,7,8),((5,9),6,7,8),((9,10),6,7,8),((10,11),6,7,8),((1,11),6,7,8)]
#letti=[((1,2),6),((2,3),6),((3,4),6),((4,5),6),((5,9),6),((9,10),6),((10,11),6),((1,11),6)]
expsi=[-1,-1,-1,-1,-1,-1,-1,1]
'''

'''
while True:
    red=reduction2(letti,expsi,gens)
    letti,expsi=red[0],red[1]
    if all([el in gens for el in letti]):
        break
print [('g'+str(gens.index(letti[i])),expsi[i]) for i in xrange(len(letti))]
print gens
'''   

#cycles for theta4
'''
#alphaU
letti=[((1,2),9,10),((2,3),9,10),((3,4),9,10),((4,5),9,10),((5,6),9,10),((1,6),9,10)]
expsi=[-1,-1,-1,-1,-1,1]
'''
'''
#alphaM
letti=[((4,5),2,3),((5,6),2,3),((1,6),2,3),((1,8),2,3),((7,8),2,3),((4,7),2,3)]
expsi=[-1,-1,1,-1,1,1]
'''
'''
#alphaD
letti=[((4,7),2,3),((7,8),2,3),((1,8),2,3),((1,10),2,3),((9,10),2,3),((4,9),2,3)]
expsi=[-1,-1,1,-1,1,1]
'''
'''
#gamma
#letti=[((4,7),1,3),((4, 5),1,3),((3, 4),1,5),((4, 7),1,5),((4, 5),1,7),((3, 4),1,7)]
#letti=[((4,9),1,3),((4, 5),1,3),((3, 4),1,5),((4, 9),1,5),((4, 5),1,9),((3, 4),1,9)]
#letti=[((4,9),1,3),((4, 7),1,3),((3, 4),1,7),((4, 9),1,7),((4, 7),1,9),((3, 4),1,9)]
letti=[((4,9),3,10),((4, 5),3,10),((3, 4),5,10),((4, 9),5,10),((4, 5),9,10),((3, 4),9,10)]
expsi=[1,-1,-1,-1,1,1]
'''
'''
while True:
    red=reduction2(letti,expsi,gens)
    letti,expsi=red[0],red[1]
    if all([el in gens for el in letti]):
        break
print [('g'+str(gens.index(letti[i])),expsi[i]) for i in xrange(len(letti))]
'''  
#cycles for K4
'''
#gamma
#letti=[((7,9),6),((7, 8),6),((6, 7),8),((7, 9),8),((7, 8),9),((6, 7),9)]
letti=[((5,10),4),((5, 6),4),((4, 5),6),((5, 10),6),((5, 6),9),((4, 5),10)]
expsi=[1,-1,-1,-1,1,1]
'''
'''
#omega1
letti=[((1,2),4),((2,3),4),((3,8),4),((7,8),4),((7,9),4),((1,9),4)]
expsi=[-1,-1,-1,1,-1,1]
'''
'''
#omega2
letti=[((3,4),10),((4,5),10),((5,6),10),((6,7),10),((7,8),10),((3,8),10)]
expsi=[-1,-1,-1,-1,-1,1]
'''
'''
#omega3
letti=[((1,9),2),((7,9),2),((6,7),2),((5,6),2),((5,10),2),((1,10),2)]
expsi=[-1,1,1,1,-1,1]
'''
'''
#omega0
letti=[((1,2),6),((2,3),6),((3,4),6),((4,5),6),((5,10),6),((1,10),6)]
expsi=[-1,-1,-1,-1,-1,1]
'''
'''
#letti=[((1,2),6,7),((2,3),6,7),((3,4),6,7),((4,5),6,7),((5,9),6,7),((9,10),6,7),((10,11),6,7),((1,11),6,7)]
letti=[((1,2),6),((2,3),6),((3,4),6),((4,5),6),((5,9),6),((9,10),6),((10,11),6),((1,11),6)]
expsi=[-1,-1,-1,-1,-1,-1,-1,1]
while True:
    red=reduction2(letti,expsi,gens)
    letti,expsi=red[0],red[1]
    if all([el in gens for el in letti]):
        break
print [('g'+str(gens.index(letti[i])),expsi[i]) for i in xrange(len(letti))]
'''