#!/usr/bin/env sage -python

import sys
from ga import GA
from sage.all import *
import numpy as np
from numpy.random import randint, rand

def fit(g):
    if g.order() < 1:
        print("empty graph")
    return g.lovasz_theta() / (len(g.independent_set()))

def mu(g):
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

def cr(g1, g2):
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



def rand_graph(n, m):
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

n = 10
pop_size = 200
pop = [rand_graph(10, randint(n, n*(n-1)/2 + 1)) for _ in range(pop_size)]

ga = GA(fit, mu, cr, 0.3, 0.1)
results = ga.run(pop, 100000, 1.15)
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
