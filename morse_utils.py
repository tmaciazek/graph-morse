import numpy as np
import os
import itertools as it
import copy
from operator import itemgetter

from graph_utils import *


def find_crit_cells(graph, dim=1, Npart=2):
    """
    Finds the critical cells of the Morse complex.

    Arguments:
    graph - adjacency matrix of the graph with deleted edges marked by '-1'
    dim - dimension of the critical cells
    Npart - number of particles

    Returns:
    cells - the set of critical cells of dimension dim of graph with Npart particles
    """
    tree, blocked_edges, del_edges, vert = group_edges(graph)
    cells = set()
    n_del_min = 2 * dim - Npart
    if n_del_min < 0:
        n_del_min = 0
    for n_del in range(n_del_min, dim + 1, 1):
        del_coll = it.combinations(del_edges, n_del)
        blocked_coll = it.combinations(blocked_edges, dim - n_del)
        for elem in it.product(del_coll, blocked_coll):
            flat = []
            for i in range(n_del):
                flat.append(elem[0][i][0])
                flat.append(elem[0][i][1])
            for i in range(dim - n_del):
                flat.append(elem[1][i][0][0])
                flat.append(elem[1][i][0][1])
                flat.append(elem[1][i][1])
            if len(set(flat)) == len(flat):
                for comb_verts in it.combinations(
                    tuple(set(vert) - set(flat)), Npart + n_del - 2 * dim
                ):
                    all_blocked = True
                    for v in comb_verts:
                        if v == 1 or (
                            next((e[0] for e in tree if e[1] == v), None)
                            in (flat + list(comb_verts))
                        ):
                            continue
                        all_blocked = False
                        break
                    if all_blocked:
                        candidate = elem[0]
                        candidate_v = comb_verts
                        for item in elem[1]:
                            candidate += (item[0],)
                            candidate_v += (item[1],)
                        cells.add(sorted_cell(candidate + candidate_v, dim))
    return cells


def boundary2(cell):
    """
    Returns the boundary word of cell.
    """
    letters = [None] * 4
    exps = [None] * 4
    letters[3] = cell[:0] + cell[1:2] + tuple(sorted(cell[2:] + (cell[0][1],)))
    exps[3] = 1
    letters[1] = cell[:0] + cell[1:2] + tuple(sorted(cell[2:] + (cell[0][0],)))
    exps[1] = -1
    letters[2] = cell[:1] + cell[2:2] + tuple(sorted(cell[2:] + (cell[1][1],)))
    exps[2] = -1
    letters[0] = cell[:1] + cell[2:2] + tuple(sorted(cell[2:] + (cell[1][0],)))
    exps[0] = 1
    return (letters, exps)


def sorted_cell(cell, dim):
    """
    Sorts elements of the cell.

    Arguments:
    cell - a cell in the discrete configuration space
    dim - dimension of the cell
    """
    return tuple(sorted(cell[:dim], key=itemgetter(0, 1))) + tuple(sorted(cell[dim:]))


def no_order_resp_edge(cell, dim, blocked_edges, del_edges):
    """
    Checks if all the edges in the given cell are non-order respecting.

    Arguments:
    cell - a cell in the discrete configuration space
    dim - dimension of the cell
    blocked_edges - list fo blocked edges
    del_edges  - list of deleted edges

    Returns:
    (True, None) - if all the edges in cell are non-order respecting
    (False, e) - if there exist an order-respecting edge in cell
    e - the minimal order respecting edge

    """
    edges_ordered = sorted(cell[:dim], key=itemgetter(1, 0))
    for ic in range(dim):
        if edges_ordered[ic] in del_edges:
            continue
        isblocked = False
        for cb in blocked_edges:
            if cb[0] == edges_ordered[ic] and cb[1] in cell[dim:]:
                isblocked = True
                break
        if isblocked == False:
            return (False, edges_ordered[ic])  # returns the minimal o.r. edge
    return (True, None)


def all_vert_blocked(cell, dim, tree):
    """
    Checks if all the vertices in the given cell are blocked.

    Arguments:
    cell - a cell in the discrete configuration space
    dim - dimension of the cell
    tree - the given spanning tree of the graph

    Returns:
    (True, None) - if all vertices in cell are blocked
    (False, v) - if there exist an unblocked vetex in cell
    v - the smallest unblocked vertex

    """
    flat = tuple([item for sublist in cell[:dim] for item in sublist]) + cell[dim:]
    for v in cell[dim:]:
        if v == 1 or next((e[0] for e in tree if e[1] == v), None) in flat:
            continue
        return (False, v)  # returns smallest unblocked vertex
    return (True, None)


def principal_reduction(cell, dim, v, tree):
    """
    Performs the principal reduction.

    Arguments:
    cell - a cell in the discrete configuration space
    dim - dimension of the cell
    v - vertex with respect to which the principal reduction is to be done

    Returns:
    reduced - the cell which is the result of the principal reduction
    """
    reduced = cell[:dim]
    reduced += (next(e for e in tree if e[1] == v),)
    reduced = sorted_cell(reduced, dim + 1)
    ind = cell.index(v)
    reduced += cell[dim:ind] + cell[ind + 1 :]
    return reduced


def reduction2(letters, exps, criticals1, graph):
    """
    Finds the Morse image of the boundary word of a word of 2-cells.

    Arguments:
    letters - list of 2-cells in the word
    exps - exponents of the letters
    criticals1 - set of critical 1-cells
    graph - adjacency matrix of the graph with deleted edges marked by '-1'

    Returns:
    letters_new - the critical cells in the Morse image of the boundary word of
                    the given word
    exps_new  - the exponends of letters in letters_new
    """
    tree, blocked_edges, del_edges, vertices = group_edges(graph)
    letters_new = []
    exps_new = []
    for i in range(len(letters)):
        if letters[i] in criticals1:
            letters_new.append(letters[i])
            exps_new.append(exps[i])
        elif no_order_resp_edge(letters[i], 1, blocked_edges, del_edges)[0]:
            unblocked = all_vert_blocked(letters[i], 1, tree)
            if unblocked[0] == False:
                match = principal_reduction(letters[i], 1, unblocked[1], tree)
                bnd = boundary2(match)
                ib = bnd[0].index(letters[i])
                if bnd[1][ib] > 0:
                    lett_toinsert = [bnd[0][(ib - j - 1) % 4] for j in range(3)]
                    exps_toinsert = [-bnd[1][(ib - j - 1) % 4] for j in range(3)]
                else:
                    lett_toinsert = [bnd[0][(ib + j + 1) % 4] for j in range(3)]
                    exps_toinsert = [bnd[1][(ib + j + 1) % 4] for j in range(3)]
                if exps[i] > 0:
                    letters_new += lett_toinsert
                    exps_new += exps_toinsert
                else:
                    lett_toinsert.reverse()
                    exps_toinsert.reverse()
                    exps_toinsert = [-sgn for sgn in exps_toinsert]
                    letters_new += lett_toinsert
                    exps_new += exps_toinsert
        else:
            e = no_order_resp_edge(letters[i], 1, blocked_edges, del_edges)[1]
            flat = (
                tuple([item for sublist in (letters[i])[:1] for item in sublist])
                + (letters[i])[1:]
            )
            for v in letters[i][1:]:
                if (
                    v != 1
                    and v < e[1]
                    and (next((e[0] for e in tree if e[1] == v), None) in flat) == False
                ):
                    match = principal_reduction(letters[i], 1, v, tree)
                    bnd = boundary2(match)
                    ib = bnd[0].index(letters[i])
                    if bnd[1][ib] > 0:
                        lett_toinsert = [bnd[0][(ib - j - 1) % 4] for j in range(3)]
                        exps_toinsert = [-bnd[1][(ib - j - 1) % 4] for j in range(3)]
                    else:
                        lett_toinsert = [bnd[0][(ib + j + 1) % 4] for j in range(3)]
                        exps_toinsert = [bnd[1][(ib + j + 1) % 4] for j in range(3)]
                    if exps[i] > 0:
                        letters_new += lett_toinsert
                        exps_new += exps_toinsert
                    else:
                        lett_toinsert.reverse()
                        exps_toinsert.reverse()
                        letters_new += lett_toinsert
                        exps_new += [-sgn for sgn in exps_toinsert]
                    break
    return (letters_new, exps_new)


def graph_braid_group(graph, n_particles=2):
    """
    Arguments:

    graph - adjacency matrix of the graph with deleted edges marked by '-1'
    n_particles - number of particles

    Returns:

    gens, rels - generators and relators of the fundamental group of the corresponding
    dicerete configuration space
    """

    cells1 = find_crit_cells(graph, dim=1, Npart=n_particles)
    gens = list(cells1)
    cells2 = find_crit_cells(graph, dim=2, Npart=n_particles)

    rels = []
    for c2 in cells2:
        bndc = boundary2(c2)
        letti, expsi = bndc[0], bndc[1]
        while True:
            red = reduction2(letti, expsi, gens, graph)
            letti, expsi = red[0], red[1]
            if all([el in gens for el in letti]):
                break
        rels.append(
            [("g" + str(gens.index(letti[i])), expsi[i]) for i in range(len(letti))]
        )

    return gens, rels

def find_morse_image(letts, exps, gens, graph):
	while True:
		letts, exps = reduction2(letts, exps, gens, graph)
		if all([el in gens for el in letts]):
			break
	return letts, exps
