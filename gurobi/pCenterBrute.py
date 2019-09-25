# Copyright 2019 Antonio Medrano
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Antonio Medrano

import numba as nb
import numpy as np
import itertools


# BRUTE FORCE METHODS FOR P = 2
def chunk2(A, n):
    global_best = np.inf

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

    return global_best, save_rows

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
def nbSerial2(A, n):
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


@nb.njit(parallel=True)
def nbParallel2(A, n):
  all_global_best=np.inf
  rows=np.empty((2),dtype=np.uint64)
  save_rows=np.empty((n,2),dtype=np.uint64)
  global_best_Temp=np.empty((n),dtype=A.dtype)
  global_best_Temp[:]=np.inf

  for i in range(n-1):
      for j in nb.prange(i+1, n):
          row_1=0
          row_2=0
          global_best=np.inf
          # find the maximum of the element-wise minimum of the three vectors
          local_best = max_min_2(A[i,:], A[j,:])
          # if local_best is lower than global_best, update global_best
          if (local_best < global_best):
              global_best = local_best
              row_1 = i
              row_2 = j

          save_rows[j,0]=row_1
          save_rows[j,1]=row_2
          global_best_Temp[j]=global_best

      ind=np.argmin(global_best_Temp)
      if (global_best_Temp[ind]<all_global_best):
          # rows=save_rows[ind,:]
          rows[0] = save_rows[ind,0]
          rows[1] = save_rows[ind,1]
          all_global_best=global_best_Temp[ind]

  return all_global_best, rows
  
  
  
# BRUTE FORCE METHODS FOR P = 3
def chunk3(A, n):
    global_best = np.inf

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

    return global_best, save_rows


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
def nbSerial3(A, n):
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


@nb.njit(parallel=True)
def nbParallel3(A, n):
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
  
  
 
# BRUTE FORCE METHODS FOR P = 4
def chunk4(A, n):
    global_best = np.inf

    chunk = 1000 # define chunk lenght, if to small, the code won't take advantatge
                 # of vectorization, if it is too large, excessive memory usage will
                 # slow down execution, or Memory Error will be risen
    combinations = itertools.combinations(range(n),4) # generate iterator containing
                                            # all possible combinations of 3 columns
    N = n*(n-1)*(n-2)*(n-3)//24 # number of combinations (length of combinations cannot be
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
        current_comb = np.fromiter(combinations,dtype='i,i,i,i',count=counts)\
                         .view(('i',4))
        chunk_best = A[current_comb].min(axis=1).max(axis=1) # maximum of element-wise
                                                             # minimum in current chunk
        ravel_save_row = chunk_best.argmin() # minimum of maximums in current chunk
        # check if current chunk contains global minimum
        if chunk_best[ravel_save_row] < global_best:
            global_best = chunk_best[ravel_save_row]
            save_rows = current_comb[ravel_save_row]
            
    return global_best, save_rows
  

@nb.njit()
def max_min_4(A,B,C,D):
  max_of_min=-np.inf
  for i in range(A.shape[0]):
    loc_min=A[i]
    if (B[i]<loc_min):
      loc_min=B[i]
    if (C[i]<loc_min):
      loc_min=C[i]
    if (D[i]<loc_min):
      loc_min=D[i]

    if (max_of_min<loc_min):
      max_of_min=loc_min

  return max_of_min


@nb.njit()
def nbSerial4(A, n):
  save_rows=np.zeros(4,dtype=np.uint64)
  global_best=np.inf
  for i in range(n-3):
      for j in range(i+1, n-2):
          for k in range(j+1, n-1):
              for l in range(k+1, n):
                  # find the maximum of the element-wise minimum of the three vectors
                  local_best = max_min_4(A[i,:], A[j,:], A[k,:], A[l,:])
                  # if local_best is lower than global_best, update global_best
                  if (local_best < global_best):
                      global_best = local_best
                      save_rows[0] = i
                      save_rows[1] = j
                      save_rows[2] = k
                      save_rows[3] = l

  return global_best, save_rows


@nb.njit(parallel=True)
def nbParallel4(A, n):
  all_global_best=np.inf
  rows=np.empty((4),dtype=np.uint64)
  save_rows=np.empty((n,4),dtype=np.uint64)
  global_best_Temp=np.empty((n),dtype=A.dtype)
  global_best_Temp[:]=np.inf

  for i in range(n-3):
      for j in nb.prange(i+1, n-2):
          global_best=np.inf
          row_1=0
          row_2=0
          row_3=0
          row_4=0
          for k in range(j+1, n-1):
              for l in range(k+1, n):
                  # find the maximum of the element-wise minimum of the three vectors
                  local_best = max_min_4(A[i,:], A[j,:], A[k,:], A[l,:])
                  # if local_best is lower than global_best, update global_best
                  if (local_best < global_best):
                      global_best = local_best
                      row_1 = i
                      row_2 = j
                      row_3 = k
                      row_4 = l

          save_rows[j,0]=row_1
          save_rows[j,1]=row_2
          save_rows[j,2]=row_3
          save_rows[j,3]=row_4
          global_best_Temp[j]=global_best

      ind=np.argmin(global_best_Temp)
      if (global_best_Temp[ind]<all_global_best):
          rows[0] = save_rows[ind,0]
          rows[1] = save_rows[ind,1]
          rows[2] = save_rows[ind,2]
          rows[3] = save_rows[ind,3]
          all_global_best=global_best_Temp[ind]

  return all_global_best, rows