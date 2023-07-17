import numpy as np
import os
import itertools as it


def load_graph_from_dataset(graph_name, subdivision):
    """
    Loads graph graph_name from the dataset in directory graphs_dataset.
    subdivison - integer specifying how 'sufficeintly subdivided'
    the graph is supposed to be

    Returns:

    Adjacency matrix of the subdivided graph with the deleted edges marked by '-1'.
    """
    available_graphs = set(os.listdir("graphs_dataset"))
    filename = graph_name + "_N" + str(subdivision)
    if filename not in available_graphs:
        print("Graph " + str(filename) + " is not in the dataset.")
        return None
    else:
        return np.loadtxt("graphs_dataset/" + filename, dtype=int, delimiter=" ")


def export_graph_to_dataset(graph_mtx, graph_name, subdivision):
    """
    Exports graph to the dataset

    graph_name - graph name
    graph_mtx - adjacency matrix of the subdivided graph
                            with the deleted edges marked by '-1'.
    subdivison - integer specifying how 'sufficeintly subdivided'
    the graph is supposed to be
    """
    filename = graph_name + "_N" + str(subdivision)
    np.savetxt("graphs_dataset/" + filename, graph_mtx, fmt="%d")
    print(
        "Graph "
        + filename.split("_N")[0]
        + " subdiv "
        + filename.split("_N")[1]
        + " has been exported."
    )

    return None


def print_graphs_from_dataset():
    print(sorted(set(os.listdir("graphs_dataset"))))

    return None


def group_edges(graph):
    """
    Groups edges of the graph for discrete Morse theory.

    Returns:

    tree - list of edges contained in the spanning tree
    blocked_edges - (technical) list of critical 1-cells in the Morse complex
                                    for N=2 particles
    del_edges - edges that do not belong to the spanning tree
    vertices - list of vertices
    """
    size = graph.shape[0]
    tree = []
    del_edges = []
    blocked_edges = []
    vertices = [i + 1 for i in range(size)]
    for com in it.combinations(range(size), 2):
        if graph[com[0]][com[1]] > 0:
            tree.append((com[0] + 1, com[1] + 1))
        elif graph[com[0]][com[1]] < 0:
            del_edges.append((com[0] + 1, com[1] + 1))
    for i in range(size):
        if (graph[i] > 0).sum() > 2:
            blocked = []
            for e in tree:
                if e[0] == vertices[i]:
                    blocked.append(e)
            for comb in it.combinations(blocked, 2):
                blocked_edges.append((comb[1], comb[0][1]))
    return tree, blocked_edges, del_edges, vertices
