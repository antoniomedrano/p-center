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

import sys
import time
import numpy as np
import itertools
import readDataFiles
import plot
from scipy.spatial.distance import cdist

def Run_pCenter(p):
    
    """ Example of solving a 2-Center problem using brute force"""  
    
    start_time = time.time()
    
    distMatrix = computeDistanceMatrix()
    #print distMatrix
    
    SDmin, locations = SolveModel(p, distMatrix)
    
    total_time = time.time()-start_time
    #SDmin = m.objVal
    
    displaySolution(locations, p, SDmin**0.5, total_time)
    
    
def computeDistanceMatrix():
        
    #declare a couple variables
    global siteIDs
    
    # Pull out just the site/demand IDs from the data
    siteIDs = sites[:,0]
    
    # Pull out just the coordinates from the data
    xyPointArray = sites[:,[1,2]]
    #A = [xyPointArray[i][:] for i in demandIDs]
    #B = [xyPointArray[j][:] for j in siteIDs]
    A = xyPointArray
    B = A
    #print A
    
    # Compute the distance matrix, using the squared distance
    distMatrix = cdist(A, B,'sqeuclidean')

    return distMatrix


def SolveModel(p, A):

    n = len(A)
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
    
    
def displaySolution(locations, p, zbest, total_time):

    print('Total problem solved in %f seconds' % total_time)
    print()
    # The objective value of the solution.
    print('p = %d' % p)
    print('SD = %f' % zbest)
    # print the selected sites
    print()
    for j in locations:
        print("Site selected %s" % int(siteIDs[j]))
    
    # plot solution
    # plot.plotSolution(sites, Y, range(numSites), SDmin)
    

def read_problem(file):
    global numSites
    global numDemands
    global sites
        
    try:
        if (file[-3:].lower() == "dat"):
            sites = readDataFiles.readDat(file)
        elif (file[-3:].lower() == "tsp"):
            sites = readDataFiles.readTSP(file)
    except IOError:
        print('Error reading file')
        raise
        
    numSites = sites.shape[0]    
    numDemands = numSites
    
    #plot.plotData(sites)
    
    print('%d locations' % numSites)
    print('Finished Reading File!')


def main(unused_argv):
    print('---- 2-Center solved via brute force -----')
    Run_pCenter(p)


""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 2 and len(sys.argv) <= 3:
    file = '../data/' + sys.argv[2]
    p = float(sys.argv[1])
    print("Problem instance from: ", file)
    read_problem(file)
    main(None)
  elif len(sys.argv) > 1 and len(sys.argv) <= 2:
    p = float(sys.argv[1])
    file = '../data/swain.dat'
    print("Problem instance from: ", file)
    read_problem(file)
    main(None)
  else:
    print("Please Pass: Service Distance; Data to Use")
    print("Problem not executed!")