import numba as nb
import numpy as np
import time

#Min and max library calls may be costly for only 3 values
@nb.njit()
def min_3(A,B,C):
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
  

@nb.njit(parallel=True)
def your_func(A):
  n=A.shape[0]
  all_global_best=np.inf
  rows=np.empty((3),dtype=np.uint64)

  save_rows=np.empty((n,3),dtype=np.uint64)
  global_best_Temp=np.empty((n),dtype=A.dtype)
  global_best_Temp[:]=np.inf

  for i in range(n):
      for j in nb.prange(i+1, n):
          global_best=np.inf
          for k in range(j+1, n):
              # find the maximum of the element-wise minimum of the three vectors
              local_best = min_3(A[i,:], A[j,:], A[k,:])
              # if local_best is lower than global_best, update global_best
              if (local_best < global_best):
                  global_best = local_best
                  row_1 = i
                  row_2 = j
                  row_3 = k

          save_rows[j,0]=row_1
          save_rows[j,1]=row_2
          save_rows[j,2]=row_3
          global_best_Temp[j]=global_best

      ind=np.argmin(global_best_Temp)
      if (global_best_Temp[ind]<all_global_best):
          rows=save_rows[ind,:]
          all_global_best=global_best_Temp[ind]

  return all_global_best, rows

  
n = 439
print "n = %d" % n
np.random.seed(2)
A = np.random.rand(n,n)
start_time = time.time()

global_best, save_rows = your_func(A)
total_time = time.time()-start_time
print global_best, save_rows
print total_time