# Copyright 2017 Antonio Medrano
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
import readDataFiles
import plot
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist
from gurobipy import *

def RunLSCP(SD):
    
    # Example of simple LSCP program with the C++ style API.
    m = Model()
    
    #print sites
    #print np.shape(sites)
    start_time = time.time()
    
    computeCoverageMatrix(SD)
    
    BuildModel(m)
    SolveModel(m)
    
    total_time = time.time()-start_time
    
    p = m.objVal
    displaySolution(m, p, total_time)
    
    
def computeCoverageMatrix(SD):
        
    #declare a couple variables
    global numDemands
    global numSites
    global cover_rows
    global cols
    global siteIDs
    
    # Pull out just the site/demand IDs from the data
    siteIDs = sites[:,0]
    
    # Pull out just the coordinates from the data
    xyPointArray = sites[:,[1,2]]
    #A = [xyPointArray[i][:] for i in demandIDs]
    #B = [xyPointArray[j][:] for j in siteIDs]
    A = xyPointArray
    #print A
    
    # Compute the distance matrix, using the squared distance
    sqDistMatrix = cdist(A, A,'sqeuclidean')
    #distances = np.unique(sqDistMatrix)
    SDsquared = SD*SD

    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)
    # Determine neighborhood of sites within 2*SD of sites (symmetric)
    C2 = (sqDistMatrix <= 4*SDsquared).astype(int)
    # NOTE: For non-symmetric problems, need to make a demand-to-demand matrix as well
    
    start_time = time.time()
    C, rows, cols = dominationTrim(C, C2)
    print 'Domination time = %f' % (time.time()-start_time)

    # shorten the facility data sets
    siteIDs = siteIDs[cols]
    numSites = len(siteIDs)
    numDemands = len(rows)

    # Convert coverage to sparse matrix
    cover_rows = [np.nonzero(t)[0] for t in C]

    return 0


def dominationTrim(A, A2):
    
    r,c = A.shape
    
    # lower triangle of coverage matrix for checking only columns within 2*SD
    # Explanation:
    # looking down each column, each row with a 1 represents a site within 2*SD of that site
    # using tril means you don't check backwards
    # NOTE: For non-symmetric problems, U should use a demand-to-demand matrix
    L = np.tril(A2,-1)
    U = np.triu(A2,1)
    
    rows = np.array(range(r))
    cols = np.array(range(c))
    
    k = 0
    
    while True:
        k += 1
        c_keeps = np.ones(c)
        r_keeps = np.ones(r)

        # create a list of sets containing the indices of non-zero elements of each column
        C = csc_matrix(A)
        D = [set(C.indices[C.indptr[i]:C.indptr[i+1]]) for i in range(len(C.indptr)-1)]

        # Column domination
        # find subsets, ignoring columns that are known to already be subsets
        for i in range(c):
            if c_keeps[i]==0:
                continue
            col1 = D[i]
            for j in np.nonzero(L[:,i])[0]:
                col2 = D[j]
                if col2.issubset(col1):
                    c_keeps[j] = 0
                elif col1.issubset(col2):
                    c_keeps[i] = 0
                    break
                
        A = A[:,c_keeps.astype(bool)]
        cols = cols[c_keeps.astype(bool)]
    
        # Row Domination
        # find subsets, ignoring rows that are known to already be subsets
        # create a list of sets containing the indices of non-zero elements of each column
        R = csr_matrix(A)
        S = [set(R.indices[R.indptr[i]:R.indptr[i+1]]) for i in range(len(R.indptr)-1)]
    
        for i in range(r):
            if r_keeps[i]==0:
                continue
            row1 = S[i]
            for j in np.nonzero(U[i,:])[0]:
                row2 = S[j]
                if row2.issubset(row1):
                    r_keeps[i] = 0
                    break
                elif row1.issubset(row2):
                    r_keeps[j] = 0
                
        A = A[r_keeps.astype(bool),:]
        rows = rows[r_keeps.astype(bool)]
        
        # Check if there was an improvement. If so, repeat.
        rnew,cnew = A.shape
        print k, rnew, cnew
                
        if (rnew == r and cnew == c):
            break
        else:
            # remaining sites to sites distance matrix
            L = L[np.ix_(c_keeps.astype(bool), c_keeps.astype(bool))]
            # remaining demands to demands distance matrix
            U = U[np.ix_(r_keeps.astype(bool), r_keeps.astype(bool))]
            r = rnew
            c = cnew            
    
    return A, rows, cols
    

def BuildModel(m):
    
    # DECLARE VARIABLES:
    # Facility Site binary decision variables X
    # Each has a coefficient of 1 in the objective
    sitesRange = range(numSites)
    X = m.addVars(sitesRange,
                  vtype=GRB.BINARY,
                  obj=np.ones(numSites),
                  name="X")
    
    # Define Coverage Constraints:
    for i in range(numDemands):
        m.addConstr(quicksum(X[j]  for  j  in  cover_rows[i])  >=  1)
    
    # The objective is to minimize the number of located facilities
    m.modelSense = GRB.MINIMIZE
    m.update()
    
    print 'Number of variables = %d' % m.numintvars
    print 'Number of constraints = %d' % m.numconstrs
    print
    return 0


def SolveModel(m):
    """Solve the problem and print the solution."""
    m.Params.OutputFlag = 0
    m.Params.ResultFile = "output.sol"
    m.optimize()
    
    
def displaySolution(m, p, total_time):

    print 'Total problem solved in %f seconds' % total_time
    print
    # The objective value of the solution.
    print 'p = %d' % p
    print 'SD = %f' % SD
    # print the selected sites
    print
    j = 0    
    for v in m.getVars():
        if (v.x == 1.0):
            print "Site selected %s" % int(siteIDs[j])
        j += 1
        
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # plot solution 
    # plot.plotSolution(sites, X, range(numSites), SD)
            
            
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
        print 'Error reading file'
        raise
        
    numSites = sites.shape[0]    
    numDemands = numSites
    
    # plot.plotData(sites)
    
    print '%d locations' % numSites
    print 'Finished Reading File!'


def RunGurobi_LSCP(SD):
    print ('---- LSCPdomdom with Gurobi -----')
    RunLSCP(SD)


def main(unused_argv):
    RunGurobi_LSCP(SD)


""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 2 and len(sys.argv) <= 3:
    file = '../data/' + sys.argv[2]
    SD = float(sys.argv[1])
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  elif len(sys.argv) > 1 and len(sys.argv) <= 2:
    SD = float(sys.argv[1])
    file = '../data/swain.dat'
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  else:
    print "Please Pass: Service Distance; Data to Use"
    print "Problem not executed!"