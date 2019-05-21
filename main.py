#!/usr/bin/env sage -python

import sys
from ga import GA
from sage.all import *
import numpy as np
import itertools
from numpy.random import randint, rand
from sage.graphs.graph_generators_pyx import RandomGNP

def fit(g):
    if g.order() < 1:
        print("empty graph")
    #return g.lovasz_theta() / (len(g.independent_set()) * g.order())
    return g.lovasz_theta() / len(g.independent_set())

def mu(g):
    """Choose a random edge uv, if exists remove it if not add it"""
    if g.order()==1:
        return RandomGNP(10, .5, directed=False)
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
        if not g.has_edge(uv):
            g.add_edge(uv)
            i+=1
    return g

def cr2(g1, g2):
    """Create a new graph by randomly sampling the product of the parents uniformly."""
    if g1.order() > 30 or g2.order() > 30:
        print "too large"
        return Graph({0:[]})
    #product = g1.strong_product(g2)
    product = g1.disjunctive_product(g2)
    #prob = 2*(g1.order()*g2.order())**-0.5
    #prob = min (1, 50/float((g1.order()*g2.order())) )
    prob = 1.0/ (len(g1.independent_set())*len(g2.independent_set()))
    #print prob
    sample = product.random_subgraph(prob)
    # while sample.order()< 7:
    #     sample = product.random_subgraph(prob)
    # print sample.order()
    if sample.order()==0:
        # print "nooo"
        # print prob
        # print product.order()
        return Graph({0:[]})
    return sample

def cr3(g1,g2):
    """Adds edges randomly between the two graphs"""
    new_graph = g1+g2
    #print new_graph.edges()
    for a,b in itertools.product(g1.vertices(),g2.vertices()):
        r = np.random.rand()
        if r < 0.5:
            new_graph.add_edge(((0,a,None),(1,b,None)))
    new_graph.random_subgraph(0.5, inplace = True)
    while new_graph.order()>50:
        new_graph.random_subgraph(0.2, inplace = True)
    if new_graph.order() ==0:
        print "too small"
        return Graph({0:[]})
    return new_graph

def cr(g1,g2):
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

n = 7 # graph size
pop_size = 100
threshold = 0.11
pop = [rand_graph(10, randint(n, n*(n-1)/2 + 1)) for _ in range(pop_size)]

ga = GA(fit, mu, cr, 0.5, 0.1)
results = ga.run(pop, 100, threshold)
results = sorted(results, key = lambda x: -x[1])
for g, fit in results:
    print(g.adjacency_matrix())
    print(g.lovasz_theta())
    print(len(g.independent_set()))
    r = g.lovasz_theta() / (len(g.independent_set()))
    print(r)
    print(fit)
    print("---------------------------------------")


# G = rand_graph(5, 6)
# print(G.edges())
# G.add_edge(1, 5)
# print(G.edges())
# # G.plot().show()
# print G.lovasz_theta()
# print len(G.independent_set())
# print G.chromatic_number()
