import numpy as np
from itertools import combinations
import time

n = 431
np.random.seed(2)
A = np.random.rand(n,n)
global_best = np.inf
#print(A)

start_time = time.time()

for i in range(len(A)-2):
    for j in range(i+1,len(A)-1):
            # find the maximum of the element-wise minimum of the three vectors
            local_best = np.amax(np.array([A[i,:], A[j,:]]).min(0))
            # if local_best is lower than global_best, update global_best
            if (local_best < global_best):
                global_best = local_best
                save_rows = np.array([i, j])

total_time = time.time()-start_time

print global_best, save_rows
print total_time

global_best = np.inf
start_time = time.time()

for i, j in combinations(range(n), 2):
    local_best = np.amax(np.array([A[i,:], A[j,:]]).min(0))
    if local_best < global_best:
        global_best = local_best
        save_rows = np.array([i, j])

total_time = time.time()-start_time

print global_best, save_rows
print total_time

global_best = np.inf
start_time = time.time()

coms = np.fromiter(combinations(np.arange(n), 2), 'i,i').view(('i', 2))
best = A[coms].min(1).max(1)
at = best.argmin()
global_best = best[at]
save_rows = coms[at]

total_time = time.time()-start_time
print global_best, save_rows
print total_time