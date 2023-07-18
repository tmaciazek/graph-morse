# Discrete Morse theory for graph configuration spaces

Python implementations of algorithms computing topological properties of unordered graph configuration spaces via the discrete Morse theory. The algorithms and the underlying theory are due to Daniel Farley and Lucas Sabalka [Algebraic & Geometric Topology Volume 5 (2005) 1075–1109, Forum Math. 24 (2012) 827-859] and Ki Hyoung Ko and Hyo Won Park  [Discrete Comput Geom (2012) 48:915–963]. 

In particular, we implement computations relating to presentations of graph braid groups.

# Folders and files in this repo

-  [tutorial.ipynb](tutorial.ipynb): Jupyter notebook showing how to use the code,
-  [morse_utils.py](morse_utils.py): utility functions relating to discrete Morse theory (implemenation of the principal reduction, finding the critical cells, implementation of boundary maps),
-  [graph_utils.py](graph_utils.py): utility functions for manipulating graphs (a graph and its spanning tree are endoced as adjacency matrix with '-1' entry whenever the corresponsing edge is a **deleted edge**, i.e. an edge which is not contained in the spanning tree).
-  [graphs_dataset](graphs_dataset): folder containing the database of graphs (the name format is `name_of_the_graph_N_subdivison`, where `subdivision` is the number of particles for which the graph is **suffciently subdivided** -- see [Abrams, A.: Configuration spaces and braid groups of graphs, PhD thesis, University of California, Berkeley, 2000]).

