#!/usr/bin/env sage -python

import sys
from ga import GA
import functions as FUN
from numpy.random import randint, rand

n = 10 # graph size
pop_size = 1000
threshold = 0.110
pop = [FUN.rand_graph(n, randint(n, n*(n-1)/2 + 1)) for _ in range(pop_size)]

ga = GA(FUN.fit_eigen_values, FUN.mu, FUN.cr4, 0.3, 0.2)
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
