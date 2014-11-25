from math import radians, cos, sin, asin, sqrt
import numpy as np

unique_test = lambda vect: len(set(np.unique([int(a) for a in vect])))

def genetic_optimization(x0,
                         cost,
                         crossover,
                         mutation,
                         generations=1e5,
                         mortality=.5,
                         mutation_probability=.1):
    """
    Generic function for maximiation via genetic algorithm.
    
    x0: Initializaiton matrix. solution dimension x population size
    cost: cost function. Takes a vector as input, returns a scalar
    crossover: crossover function. Takes two vectors as inputs (parents), 
        returns a vector that is a convolution of the inputs (child)
    mutation: takes a vector as input, returns a vector that is a modification of the input.
    population: number of concurrents "organisms."
    generations: number of iterations to run the algorithm.
    dim: dimensions of the solution vector.
    """
    n,m = x0.shape # n = solution dimension, m = population size
    k=np.ceil(m*mortality)
    solution = {'value':x0[:,0], 'fitness':cost(x0[:,0])} # initialize solution with a random valid entry
    for j in range(generations):
        fitness = np.apply_along_axis(cost, 0, x0)        
        #survivors = fitness.argsort()[:k:-1]
        survivors = fitness.argsort()[-k:]
        if fitness[survivors[-1]] > solution['fitness']:
            solution = {'value':x0[:,survivors[-1]], 'fitness':fitness[survivors[-1]]}
            print j, solution['fitness']
        fitness = fitness[survivors] # kill of 50% least fit.        
        prob = fitness/fitness.sum() # normalize (probably not necessary)
        
        #children = np.empty_like(x0) # maybe use zeros here instead of empty_like?
        children = np.zeros(x0.shape) 
        for i in range(m):
            a,b = np.random.choice(survivors, 2, replace=False, p=prob)
            v = crossover(x0[:,a], x0[:,b])
            
            if unique_test(v) <10:
                print "crossover"
                print j,i, a, b, unique_test(v)
                print v
                raise Exception
            if np.random.random() > mutation_probability:
                v = mutation(v)
                if unique_test(v) <10:
                    print "mutation"
                    print j,i, a, b, unique_test(v)
                    print v
                    raise Exception
            children[:,i] = v
        #x0 = children.copy() # Doesn't solve the problem
        x0 = children
        
    return solution
        

def propose_path(path, n=None, swaps=1, closed=True):    
    """
    Proposal function for stochastic solutions to traveling salesman. Gives a "neighbor" path
    to the in put path.
    """
    swaps = np.random.choice(range(1,swaps+1))
    if not n:
        n = len(path)
    path2 = path.copy() # we don't want to operate on the input
    if closed:
        ix = np.random.choice(range(1,n), 2*swaps) # ignore first position in closed circuit, enforce all candidates start/end at same location.
    else:
        ix = np.random.choice(range(n), 2*swaps)
    path2[ix] = path2[ix[::-1]]
    return path2 # do we need the return value here? Are we modifying in place?
    
def convolve_paths(path1, path2, n, closed=True):
    """
    Combines two "parent" paths into a "child" path by taking a random 
    continuous sequence from the first path, injecting that into the child, and
    then back filling missing path nodes into the child in the order in which 
    they appear in the second parent.
    """
    # NB: this approach is *biased*. Interval length between two discrete uniforms of length n is n/3.
    #ix = np.random.choice(path1, 2, replace=False)
    if closed:
        ix = np.random.choice(np.arange(1,n), 2, replace=False)
        child = np.zeros(n+1) # This is also a necessary part of the fix
    else:
        ix = np.random.choice(np.arange(n), 2, replace=False) # THIS fixed the problem
        child = np.zeros(n) # This is also a necessary part of the fix
    ix.sort()
    s = slice(*ix)
    #child = np.empty_like(path1)     
    
    child[s] = path1[s]
    ix2 = np.setdiff1d(np.arange(n), np.arange(*ix))
    if closed:
        i=1
    else:
        i=0
    for v in path2:
        if v not in child:
            child[ix2[i]] = v
            i+=1
            if i == len(ix2):
                break
    return child


def haversine(latlon1, latlon2):
    # http://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points    
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """

    latlon1 = np.radians(latlon1)
    latlon2 = np.radians(latlon2)
    
    dlat, dlon = latlon2 - latlon1
    lat1, lon1 = latlon1
    lat2, lon2 = latlon2

    # haversine formula 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 

    # 6367 km is the radius of the Earth
    km = 6367 * c
    return km 
    
def euclidean(v1, v2):
    return np.sqrt(((v1-v2)**2).sum())
    
def path_dist(path, places, distance_func, closed=False):
    """
    Returns negative of total distance of path. Returns negative because 
    the genetic optimization is currently parameterized as a maximization.
    
    Assumes places is a dataframe.
    """
    if closed:
        path = np.hstack([path, path[0]])
    trips = np.vstack([path[:-1], path[1:]]).T
    # There's probably a vectorized way to do this...
    res = 0
    for t in trips:
        i,j = t
        res += distance_func(places[i,:], places[j,:])
    return res
    
def initialize_path_population(pathlen, pop_size, closed=True):
    if closed:
        mat = np.empty([pathlen+1, pop_size])
    else:
        mat = np.empty([pathlen, pop_size])
    #mat = np.empty([pathlen, pop_size], dtype=int)
    #mat = np.zeros([pathlen, pop_size], dtype=int)
    for i in range(pop_size):
        if closed:
            v = np.random.choice(range(1, pathlen), pathlen-1, replace=False)
            v = np.hstack([0,v,0])
        else:
            v = np.random.choice(range(pathlen), pathlen, replace=False)
        mat[:,i] = v
    return mat
    
def genetic_salesman(places, pop_size, distance_func, swaps=1, closed=True, **kwargs):
    n=places.shape[0]
    path0 = initialize_path_population(n, pop_size, closed=closed)
    return genetic_optimization(x0 = path0,
                                cost = lambda p: -1*path_dist(p, places, distance_func, closed=False),
                                crossover = lambda x,y: convolve_paths(x,y, n, closed=closed),
                                mutation = lambda p: propose_path(p, n=n, swaps=swaps, closed=closed),
                                **kwargs
                                )

#################################################

if __name__ == '__main__':
    # Test with simple euclidean coordinates.

    import matplotlib.pyplot as plt
    np.random.seed(1)
    # generate random points and plot
    pathlen = 20
    nodes = np.random.choice(range(100), 2*pathlen, replace=True).reshape([pathlen, 2])

    #test = genetic_salesman(nodes, 500, euclidean, generations=50000, mortality=.1, swaps=1, mutation_probability=.05)
    test = genetic_salesman(nodes, 200, euclidean, generations=1000, mortality=.05, swaps=1, mutation_probability=.05)
    solution = np.asarray(test['value'], dtype=int)
    xy = nodes[solution]
    plt.scatter(*zip(*nodes))
    plt.plot(*xy.T)
    plt.show()