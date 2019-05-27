"""A library of fitness, mutation, and crossover functions."""
import sys
import numpy as np
import itertools
import bronkerbosch as BON
import random
from numpy.random import randint, rand
from sage.all import *
from sage.graphs.graph_generators_pyx import RandomGNP

"""Fitness Functions"""

def fit(g):
    if g.order() < 1:
        print("empty graph")
    return g.lovasz_theta() / len(g.independent_set())

def fit_eigen_values(g):
    """Returns the ratio between the largest and second largest abs. value eigenvectors."""
    """This doesn't give good results, because we usually must assume the graphs are regular."""
    adjacency = np.array(g.adjacency_matrix())
    eigenvalues = np.linalg.eigh(adjacency)[0]
    largest = eigenvalues[-1]
    second_largest = max(abs(eigenvalues[0]),abs(eigenvalues[-2]))
    return (largest - second_largest) / largest
#TODO: reward regularity by calculating the standard deviation of the degree list.

"""Mutation Functions"""
def mu(g):
    """Choose a random edge uv, if exists remove it if not add it"""
    g = g.copy()
    v = randint(0, g.order())
    u = randint(0, g.order())
    while u == v:
        u = randint(0, g.order())
    if g.has_edge(u, v):
        if g.size() > 1:
            g.delete_edge(u, v)
    else:
        g.add_edge(u, v)

    return g

def update_independent_sets(indep_sets, edge, neighbors0, neighbors1):
    """Assumes that edge is removed and updates the maximal independent sets."""
    new_indep_sets = []
    added_new_edge = False
    for i in indep_sets:
        if (not edge[0] in i) and (not edge[1] in i):
            new_indep_sets.append(i)
        elif edge[0] in i:
            if neighbors1.intersection(i)==set([]):
                new_indep_sets.append(i.union({edge[1]}))
                added_new_edge = True
            else:
                new_indep_sets.append(i)
        elif edge[1] in i:
            if neighbors0.intersection(i)==set([]):
                new_indep_sets.append(i.union({edge[0]}))
                added_new_edge = True
            else:
                new_indep_sets.append(i)
    if added_new_edge == False:
        new_indep_sets.append({edge[0],edge[1]})
    unique_indep_sets = []
    for i in new_indep_sets: #make unique
        if not i in unique_indep_sets:
            unique_indep_sets.append(i)
    return unique_indep_sets

def remove_extra_edges(g):
    new_graph = g.copy()
    edges = len(new_graph.edges())
    indep_sets = None
    new_graph, indep_sets = _remove_extra_edge(new_graph, indep_sets)
    while(len(new_graph.edges()) != edges ):
        edges = len(new_graph.edges())
        new_graph, indep_sets = _remove_extra_edge(new_graph, indep_sets)
    return new_graph, indep_sets

def _remove_extra_edge(g, indep_sets = None):
    """Calculates the simplical complex of independent sets of g.
    If an edge doesnt intersect a maximal independent set, it can be removed
    without increasing the size of the independence number.
    """
    dict = BON.dict_from_adjacency_matrix(g.complement())
    if indep_sets is None:
        indep_sets = BON.find_cliques(dict)
    max_size = 0
    max_indep_sets = []
    new_graph = g.copy()
    for i in indep_sets:
        if len(i) > max_size:
            max_size = len(i)
            max_indep_sets = [i]
        elif len(i) == max_size:
            max_indep_sets.append(i)
    vertices_in_max_indep_set = reduce(lambda x,y: union(x,y), max_indep_sets, set([]))
    removeable_edges = [e for e in g.edges() if (not e[0] in vertices_in_max_indep_set) and (not e[1] in vertices_in_max_indep_set)]
    if len(removeable_edges)==0:
        print "no edges to remove"
        return new_graph, indep_sets
    else:
        e = removeable_edges[randint(0,len(removeable_edges)-1)]
        print "deleting ", e
        new_graph.delete_edge(e)
        #In the future, use update independent sets instead
        new_indep_sets = BON.find_cliques((BON.dict_from_adjacency_matrix(new_graph.complement())))
        return new_graph, new_indep_sets


"""Crossover Functions"""
def cr1(g1, g2):
    """Create a new graph and add edges randomly from parents."""
    e1 = g1.edges()
    e2 = g2.edges()
    g = Graph({v:[] for v in range(0, g1.order())})
    m = (g1.size() + g2.size()) // 2

    i = 0
    while i < m:
        if rand() < 0.5:
            e = e1
        else:
            e = e2
        uv = e[randint(0, len(e))]
        g.add_edge(uv)
        i+=1
    return g

def cr2(g1, g2):
    """Create a new graph by randomly sampling the product of the parents uniformly."""
    #if not g.has_edge(uv):
    if g1.order() > 30 or g2.order() > 30:
        print "too large"
        return Graph({0:[]})
    product = g1.disjunctive_product(g2)
    prob = 1.0/ (len(g1.independent_set())*len(g2.independent_set()))
    sample = product.random_subgraph(prob)
    if sample.order()==0:
        return Graph({0:[]})
    return sample

def cr3(g1,g2,downsample = False):
    """Adds edges randomly between the disjoint union of the two graphs"""
    new_graph = g1.disjoint_union(g2, labels='pairs')
    print new_graph.vertices()
    for a,b in itertools.product(g1.vertices(),g2.vertices()):
        r = np.random.rand()
        if r < 0.5:
            new_graph.add_edge(((0,a),(1,b)))
    if downsample:
        new_graph.random_subgraph(0.5, inplace = True)
    while new_graph.order()>50:
        new_graph.random_subgraph(0.2, inplace = True)
    if new_graph.order() ==0:
        print "too small"
        return Graph({0:[]})
    return new_graph

def cr4(g1,g2):
    """Keeps edges that are in both, flips a coin for edges that are in one but not the other."""
    new_graph = g1.copy()
    for edge in set(g1.edges()) ^ set(g2.edges()) :
        r = np.random.rand()
        if r < 0.5:
            if new_graph.has_edge(edge):
                    new_graph.delete_edge(edge)
            else:
                new_graph.add_edge(edge)
    return new_graph


def rand_graph(n, m):
    "Generate a random graph with n vertices and m edges"
    g = { v: [] for v in range(n)}
    i = 0
    while i < m:
        x = randint(0, n)
        y = randint(0, n)
        if x > y:
            x, y = y, x
        if x != y and y not in g[x]:
            g[x].append(y)
            i += 1
    return Graph(g)
