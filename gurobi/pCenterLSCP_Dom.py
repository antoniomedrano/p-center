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

import time
import numpy as np
import readDataFiles
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist
from gurobipy import *
setParam('OutputFlag', 0)   # mute solver meta-info

def Run_pCenterLSCP():
    
    """Example of complete p-Center program with the Gurobi API"""
    m = Model()
    
    start_time = time.time()
    
    sqDistances, sqDistMatrix = computeDistances()
    
    # p = numSites, SD = 0 is a trivial solution
    print('  p, SD')
    p = numSites
    SDsquared = 0   
    displaySolution(p, SDsquared)
    
    solution = np.empty([numSites, 2])
    solution[:,0] = range(1, numSites+1)
    solution[p-1,1] = 0
    currP = numSites
    
    SDsquared = sqDistances[0]
    essential = computeCoverageMatrix(sqDistMatrix, SDsquared)

    BuildModel(m)
    SolveModel(m)

    # get the solution and clear the solver
    p = m.objVal
    
    # check the output
    while (p < currP):
        currP -= 1
        solution[currP-1,1] = SDsquared**0.5
        displaySolution(currP, SDsquared)
    
    for k in range(1,len(sqDistances)):
        SDsquared = sqDistances[k]
        essential = computeCoverageMatrix(sqDistMatrix, SDsquared)

        UpdateModel(m)
        SolveModel(m)

        # get the solution and clear the solver
        p = m.objVal
        
        # check the output
        while (p < currP):
            currP -= 1
            solution[currP-1,1] = SDsquared**0.5
            displaySolution(currP, SDsquared)

        # terminate the search when p == 1
        if (p == 2):
            p = 1
            SDsquared = np.amin(np.amax(sqDistMatrix,0))
            solution[p-1,1] = SDsquared**0.5
            displaySolution(p, SDsquared)
            iters = k+1
            break
        if (p == 1):
            iters = k
            break
        
    total_time = time.time()-start_time
    #print solution
    print()
    print('%d LSCP distances evaluated' % iters)
    print('Total problem solved in %f seconds' % total_time)
    print()
    #plot.plotTradeoff(solution)
    

def computeDistances():
        
    #declare a couple variables
    global distances
    global siteIDsAll
    
    # Pull out just the site/demand IDs from the data
    siteIDsAll = sites[:,0]
    
    # Pull out just the coordinates from the data
    xyPointArray = sites[:,[1,2]]

    A = xyPointArray
    B = A
    
    # Compute the distance matrix, using the squared distance
    sqDistMatrix = cdist(A, B,'sqeuclidean')
    sqDistances = np.unique(sqDistMatrix)
    
    return sqDistances, sqDistMatrix
        
    
def computeCoverageMatrix(sqDistMatrix, SDsquared):
        
    #declare a couple variables
    global numDemands
    global numSites
    global cover_rows
    global cols
    global rows
    global siteIDs

    # Pull out just the site/demand IDs from the data
    siteIDs = sites[:,0]

    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)
    # Determine neighborhood of sites within 2*SD of sites (symmetric)
    C2 = (sqDistMatrix <= 4*SDsquared).astype(int)

    # Perform column domination
    C, cols = dominationTrim(C, C2)

    #print siteIDs, cols
    # shorten the facility data sets
    siteIDs = [siteIDs[j] for j in cols]
    numSites = len(siteIDs)

    # Convert coverage to sparse matrix
    cover_rows = [np.nonzero(t)[0] for t in C]

    return 0


def dominationTrim(A, A2):
    
    r,c = A.shape
    c_keeps = np.ones(c)
    cols = np.array(range(c))
    
    # lower triangle of coverage matrix for checking only columns within 2*SD
    # Explanation:
    # looking down each column, each row with a 1 represents a site within 2*SD of that site
    # using tril means you don't check backwards
    B = np.tril(A2,-1)
    
    # start_time = time.time()
    # create a list of sets containing the indices of non-zero elements of each column
    C = csc_matrix(A)
    D = [set(C.indices[C.indptr[i]:C.indptr[i+1]]) for i in range(len(C.indptr)-1)]
    # print 'Matrix to List of Sets CSC Time = %f' % (time.time()-start_time)
    
    # find subsets, ignoring columns that are known to already be subsets
    for i in cols:
        if c_keeps[i]==0:
            continue
        col1 = D[i]
        for j in np.nonzero(B[:,i])[0]:
            # I tried `if keeps[j]==false: continue` here, but that was slower
            # if keeps[j]==False: continue
            col2 = D[j]
            if col2.issubset(col1):
                c_keeps[j] = 0
            elif col1.issubset(col2):
                c_keeps[i] = 0
                break
    
    #Z = A[np.ix_(c_keeps.astype(bool),c_keeps.astype(bool))]
    A = A[:,c_keeps.astype(bool)]
    cols = cols[c_keeps.astype(bool)]
    
    return A, cols
    

def BuildModel(m):
    
    # DECLARE VARIABLES:
    # Facility Site binary decision variables X
    # Each has a coefficient of 1 in the objective
    X = m.addVars(numSites,
                  vtype=GRB.BINARY,
                  obj=1)
    
    # Define Coverage Constraints:
    for i in range(numDemands):
        m.addConstr(quicksum(X[j] for j in cover_rows[i]) >= 1)
    
    # The objective is to minimize the number of located facilities
    m.modelSense = GRB.MINIMIZE
    
    # m.update()
    # print 'Number of variables = %d' % m.numintvars
    # print 'Number of constraints = %d' % m.numconstrs
    # print
    return 0


def UpdateModel(m):
    m.remove(m.getVars())
    m.remove(m.getConstrs())
    
    # DECLARE VARIABLES:
    # Facility Site binary decision variables X
    # Each has a coefficient of 1 in the objective
    X = m.addVars(numSites,
                  vtype=GRB.BINARY,
                  obj=np.ones(numSites),
                  name="X")

    # Define Coverage Constraints:
    for i in range(numDemands):
        m.addConstr(quicksum(X[j] for j in cover_rows[i]) >= 1)


def SolveModel(m):
    """Solve the problem and print the solution."""
    # m.Params.ResultFile = "output.sol"
    m.optimize()
    
    
def displaySolution(p, SDsquared):
    # The objective value and the minimum service distance
    print('%3d, %f' % (p, SDsquared**0.5))
    
            
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
    # plot.plotData(sites)    
    print('%d locations' % numSites)


def main(unused_argv):
    print ('---- pCenterLSCP with Gurobi -----')
    Run_pCenterLSCP()


""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 1 and len(sys.argv) <= 2:
    file = '../data/' + sys.argv[1]
    print()
    print("Problem instance from: ", file)
    read_problem(file)
    main(None)
  elif len(sys.argv) > 0 and len(sys.argv) <= 1:
    file = '../data/swain.dat'
    print()
    print("Problem instance from: ", file)
    read_problem(file)
    main(None)
  else:
    print("Please Pass: Service Distance; Data to Use")
    print("Problem not executed!")