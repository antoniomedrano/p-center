import sys
import numba as nb
import numpy as np
import itertools
import time

def main(input):
  n = int(input)
  print("n = %d" % n)
  np.random.seed(2)
  A = np.random.rand(n,n)
  #print(A)

  # global_best = np.inf
  # start_time = time.time()
  
  # for i, j, k in itertools.combinations(range(n), 3):
  #     # local_best = np.amax(np.array([A[i,:], A[j,:], A[k,:]]).min(0))
  #     # local_best = np.amax(np.minimum(np.minimum(A[i,:], A[j,:]), A[k,:]))
  #     local_best = np.minimum(np.minimum(A[i,:], A[j,:]), A[k,:]).max(0)
  #     if local_best < global_best:
  #         global_best = local_best
  #         save_rows = np.array([i, j, k])
  
  # total_time = time.time()-start_time
  
  # print('standard')
  # print(global_best, save_rows)
  # print(total_time)
  
  
  # global_best = np.inf
  # start_time = time.time()
  
  # coms = np.fromiter(itertools.combinations(np.arange(n), 3), 'i,i,i').view(('i', 3))
  # #print len(coms)
  # best = A[coms].min(1).max(1)
  # at = best.argmin()
  # global_best = best[at]
  # save_rows = coms[at]
  
  # total_time = time.time()-start_time
  
  # print('vectored')
  # print(global_best, save_rows)
  # print(total_time)


  global_best = np.inf
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

  print('chunked vectored')
  print(global_best, save_rows)
  print(total_time)


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
    for i in range(n-2):
        for j in range(i+1, n-1):
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

  start_time = time.time()
  global_best, save_rows = your_func(A)
  total_time = time.time()-start_time

  print('numba serial')
  print(global_best, save_rows)
  print(total_time)


  @nb.njit(parallel=True)
  def your_func_parallel(A):
    n=A.shape[0]
    all_global_best=np.inf
    rows=np.empty((3),dtype=np.uint64)
    save_rows=np.empty((n,3),dtype=np.uint64)
    global_best_Temp=np.empty((n),dtype=A.dtype)
    global_best_Temp[:]=np.inf

    for i in range(n-2):
        for j in nb.prange(i+1, n-1):
            row_1=0
            row_2=0
            row_3=0
            global_best=np.inf
            for k in range(j+1, n):
                # find the maximum of the element-wise minimum of the three vectors
                local_best = max_min_3(A[i,:], A[j,:], A[k,:])
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
            # rows=save_rows[ind,:]
            rows[0] = save_rows[ind,0]
            rows[1] = save_rows[ind,1]
            rows[2] = save_rows[ind,2]
            all_global_best=global_best_Temp[ind]

    return all_global_best, rows

  start_time = time.time()
  global_best, save_rows = your_func_parallel(A)
  total_time = time.time()-start_time

  print('numba parallel')
  print(global_best, save_rows)
  print(total_time)

""" Main will take in 1 arguments: n-Facilities"""
if __name__ == '__main__':
  if len(sys.argv) > 0:
    n = sys.argv[1]
    main(n)
  else:
    print("Please tell us how big of a problem to solve")
    print("Problem not executed!")