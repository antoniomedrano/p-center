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

  # for i, j in itertools.combinations(range(n), 2):
  #     local_best = np.amax(np.array([A[i,:], A[j,:]]).min(0))
  #     if local_best < global_best:
  #         global_best = local_best
  #         save_rows = np.array([i, j])

  # total_time = time.time()-start_time

  # print('standard')
  # print(global_best, save_rows)
  # print(total_time)


  # global_best = np.inf
  # start_time = time.time()

  # coms = np.fromiter(itertools.combinations(np.arange(n), 2), 'i,i').view(('i', 2))
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
  combinations = itertools.combinations(range(n),2) # generate iterator containing
                                          # all possible combinations of 3 columns
  N = (n*n-n)//2 # number of combinations (length of combinations cannot be
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
      current_comb = np.fromiter(combinations,dtype='i,i',count=counts)\
                      .view(('i',2))
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
  def max_min_2(A,B):
    max_of_min=-np.inf
    for i in range(A.shape[0]):
      loc_min=A[i]
      if (B[i]<loc_min):
        loc_min=B[i]

      if (max_of_min<loc_min):
        max_of_min=loc_min

    return max_of_min

  @nb.njit()
  def your_func(A):
    n=A.shape[0]
    save_rows=np.zeros(2,dtype=np.uint64)
    global_best=np.inf
    for i in range(n-1):
        for j in range(i+1, n):
                # find the maximum of the element-wise minimum of the three vectors
                local_best = max_min_2(A[i,:], A[j,:])
                # if local_best is lower than global_best, update global_best
                if (local_best < global_best):
                    global_best = local_best
                    save_rows[0] = i
                    save_rows[1] = j

    return global_best, save_rows
    

  start_time = time.time()
  global_best, save_rows = your_func(A)
  total_time = time.time()-start_time

  print('numba serial')
  print(global_best, save_rows)
  print(total_time)


  # @nb.njit(parallel=True)
  # def your_func_parallel(A):
  #   n=A.shape[0]
  #   all_global_best=np.inf
  #   rows=np.empty((2),dtype=np.uint64)
  #   save_rows=np.empty((n,2),dtype=np.uint64)
  #   global_best_Temp=np.empty((n),dtype=A.dtype)
  #   global_best_Temp[:]=np.inf

  #   for i in range(n-1):
  #       for j in nb.prange(i+1, n):
  #           global_best=np.inf
  #           # find the maximum of the element-wise minimum of the three vectors
  #           local_best = max_min_2(A[i,:], A[j,:])
  #           # if local_best is lower than global_best, update global_best
  #           if (local_best < global_best):
  #               global_best = local_best
  #               row_1 = i
  #               row_2 = j

  #           save_rows[j,0]=row_1
  #           save_rows[j,1]=row_2
  #           global_best_Temp[j]=global_best

  #       ind=np.argmin(global_best_Temp)
  #       if (global_best_Temp[ind]<all_global_best):
  #           # rows=save_rows[ind,:]
  #           rows[0] = save_rows[ind,0]
  #           rows[1] = save_rows[ind,1]
  #           all_global_best=global_best_Temp[ind]

  #   return all_global_best, rows

  # start_time = time.time()
  # global_best, save_rows = your_func_parallel(A)
  # total_time = time.time()-start_time

  # print 'numba parallel'
  # print global_best, save_rows
  # print total_time

""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 0:
    n = sys.argv[1]
    print(n)
    main(n)
  else:
    print("Please tell us how big of a problem to solve")
    print("Problem not executed!")