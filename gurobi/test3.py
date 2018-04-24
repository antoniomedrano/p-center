import numpy as np
import itertools
import time

n = 200
np.random.seed(2)
A = np.random.rand(n,n)
global_best = 1000000000000000.0
#print(A)

# start_time = time.time()
#
# for i in range(len(A)-2):
#     for j in range(i+1,len(A)-1):
#         for k in range(j+1,len(A)):
#             # find the maximum of the element-wise minimum of the three vectors
#             local_best = np.amax(np.array([A[i,:], A[j,:], A[k,:]]).min(0))
#             # if local_best is lower than global_best, update global_best
#             if (local_best < global_best):
#                 global_best = local_best
#                 save_rows = [i, j, k]
#
# total_time = time.time()-start_time
#
# print global_best, save_rows
# print total_time
#
#
start_time = time.time()

for i, j, k in itertools.combinations(range(n), 3):
    # local_best = np.amax(np.array([A[i,:], A[j,:], A[k,:]]).min(0))
    # local_best = np.amax(np.minimum(np.minimum(A[i,:], A[j,:]), A[k,:]))
    local_best = np.minimum(np.minimum(A[i,:], A[j,:]), A[k,:]).max(0)
    if local_best < global_best:
        global_best = local_best
        save_rows = np.array([i, j, k])

total_time = time.time()-start_time

print global_best, save_rows
print total_time


start_time = time.time()

coms = np.fromiter(itertools.combinations(np.arange(n), 3), 'i,i,i').view(('i', 3))
best = A[coms].min(1).max(1)
at = best.argmin()
global_best = best[at]
save_rows = coms[at]

total_time = time.time()-start_time
print global_best, save_rows
print total_time


start_time = time.time()

chunk = 1000 # define chunk lenght, if to small, the code won't take advantatge 
             # of vectorization, if it is too large, excessive memory usage will 
             # slow down execution, or Memory Error will be risen 
combinations = itertools.combinations(range(n),3) # generate iterator containing 
                                        # all possible combinations of 3 columns
N = n*(n-1)*(n-2)//6 # number of combinations (length of combinations cannot be 
                     # retrieved because it is an iterator)
# generate a list containing how many elements of combinations will be retrieved 
# per iteration
n_chunks, remainder = divmod(N,chunk)
counts_list = [chunk for _ in range(n_chunks)]
if remainder:
    counts_list.append(remainder)

# Iterate one chunk at a time, using vectorized code to treat the chunk
for counts in counts_list:
    # retrieve combinations in current chunk
    current_comb = np.fromiter(combinations,dtype='i,i,i',count=counts)\
                     .view(('i',3)) 
    chunk_best = A[current_comb].min(axis=1).max(axis=1) # maximum of element-wise 
                                                         # minimum in current chunk
    ravel_save_row = chunk_best.argmin() # minimum of maximums in current chunk
    # check if current chunk contains global minimum
    if chunk_best[ravel_save_row] < global_best: 
        global_best = chunk_best[ravel_save_row]
        save_rows = current_comb[ravel_save_row]

total_time = time.time()-start_time
print global_best, save_rows
print total_time