import numba as nb
import numpy as np
import time

#Min and max library calls may be costly for only 3 values
@nb.njit()
def max_min_3(A,B,C):
  max_of_min=-np.inf
  for i in range(A.shape[0]):
    loc_min=A[i]
    if (B[i]<loc_min):
      loc_min=B[i]
    if (C[i]<loc_min):
      loc_min=C[i]

    if (max_of_min<loc_min):
      max_of_min=loc_min

  return max_of_min

@nb.njit()
def your_func(A):
  n=A.shape[0]
  save_rows=np.zeros(3,dtype=np.uint64)
  global_best=np.inf
  for i in range(n):
      for j in range(i+1, n):
          for k in range(j+1, n):
              # find the maximum of the element-wise minimum of the three vectors
              local_best = max_min_3(A[i,:], A[j,:], A[k,:])
              # if local_best is lower than global_best, update global_best
              if (local_best < global_best):
                  global_best = local_best
                  save_rows[0] = i
                  save_rows[1] = j
                  save_rows[2] = k

  return global_best, save_rows

  
n = 439
print "n = %d" % n
np.random.seed(2)
A = np.random.rand(n,n)
start_time = time.time()

global_best, save_rows = your_func(A)
total_time = time.time()-start_time
print global_best, save_rows
print total_time