from ga import GA
import numpy as np

def fit(x):
    return 10 - (x[0] + x[1]) ** 2

def mu(x):
    new_x = x.copy()
    i = np.random.randint(len(x))
    new_x[i] = np.random.randint(-5, 5)
    return new_x

def cr(x1, x2):
    new_x = x1.copy()
    new_x[1] = x2[1]
    return new_x

g = GA(fit, mu, cr, 0.5, 0.1)
pop = [
        [5, 1],
        [3, 7],
        [12, 2],
        [1, 2],
        [13, 3],
        [1, 25],
        [-1, 7],
        [32, 34],
        [60, 50],
        [100, 1],
        [70, 21]
        ]
r = g.run(pop, 10000, 8)
print(r)
