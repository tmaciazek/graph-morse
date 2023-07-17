import numpy as np
import os

from graph_utils import *
from morse_utils import *

tree, blocked_edges, del_edges, vertices = group_edges(graph)

gens, rels = graph_braid_group(graph)

print("Generators:")
for i in range(len(gens)):
	print(
		'g'+str(i)+' = '+str(gens[i])
	)
print("Relators:")
for rel in rels:
	print( rel )
