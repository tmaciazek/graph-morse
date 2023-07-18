import numpy as np
import os
import matplotlib.pyplot as plt

from graph_utils import *
from morse_utils import *

import networkx as nx

graph = load_graph_from_dataset( 'K3,3', 2 )
print(graph)

tree, blocked_edges, del_edges, vertices = group_edges(graph)

name='Theta4_CMP'
N0=15
graph=np.zeros((N0,N0))
graph[0][1]=1
graph[0][14]=-1
graph[1][2]=1
graph[2][3]=1
graph[2][5]=1
graph[2][7]=1
graph[3][4]=1
graph[4][12]=-1
graph[5][6]=1
graph[6][10]=-1
graph[7][8]=1
graph[8][9]=1
graph[8][11]=1
graph[8][13]=1
graph[9][10]=1
graph[11][12]=1
graph[13][14]=1

export_graph_to_dataset(graph, name, 3)

plot_graph(graph)

plt.show()
		

