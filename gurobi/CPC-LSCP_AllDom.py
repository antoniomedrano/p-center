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
    p = m.objVal + len(essential)
    
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
        p = m.objVal + len(essential)
        
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

    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)
    # Determine neighborhood of sites within 2*SD of sites (symmetric)
    C2 = (sqDistMatrix <= 4*SDsquared).astype(int)
    # NOTE: For non-symmetric problems, need to make a demand-to-demand matrix as well

    start_time = time.time()
    C, rows, cols, essential = dominationTrim(C, C2)
    #print 'Domination time = %f' % (time.time()-start_time)

    # shorten the facility data sets
    siteIDs = siteIDsAll[cols]
    numSites = len(siteIDs)
    numDemands = len(rows)

    # Convert coverage to sparse matrix
    cover_rows = [np.nonzero(t)[0] for t in C]

    return np.nonzero(essential)[0]


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
    essential = np.zeros(c)
        
    while True:
        c_keeps = np.ones(c)
        r_keeps = np.ones(r)

        # create a list of sets containing the indices of non-zero elements of each column
        C = csc_matrix(A)
        D = [set(C.indices[C.indptr[i]:C.indptr[i+1]]) for i in range(len(C.indptr)-1)]

        # COLUMN DOMINATION
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
        # remaining sites to sites distance matrix
        L = L[np.ix_(c_keeps.astype(bool), c_keeps.astype(bool))]

        # ROW DOMINATION
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
        # remaining demands to demands distance matrix
        U = U[np.ix_(r_keeps.astype(bool), r_keeps.astype(bool))]        


        # ESSENTIAL SITES
        # Designate sites that uniquely cover certain demands as essential, requiring forced
        # location there.
        r_keeps = np.ones(len(rows))
        c_keeps = np.ones(len(cols))
        rSum = np.sum(A, axis=1)

        for i in range(len(rSum)):
            if rSum[i] == 1:
                r_keeps[i] = 0
                c_keeps[np.nonzero(A[i,:])] = 0
                essential[cols[np.nonzero(A[i,:])]] = 1
        
        A = A[np.ix_(r_keeps.astype(bool), c_keeps.astype(bool))]
        cols = cols[c_keeps.astype(bool)]
        rows = rows[r_keeps.astype(bool)]
        
        
        # CHECK IF SHOULD REPEAT
        # Check if there was an improvement. If so, repeat.
        rnew,cnew = A.shape
        #print k, rnew, cnew
        
        if (rnew == 0):
            cols = rows  # make sure no problem gets formulated and solved
            break
        if (cnew == 0):
            rows = cols  # make sure no problem gets formulated and solved
            break
        if (rnew == r and cnew == c):
            break
        else:
            L = L[np.ix_(c_keeps.astype(bool), c_keeps.astype(bool))]
            U = U[np.ix_(r_keeps.astype(bool), r_keeps.astype(bool))]        
            r = rnew
            c = cnew
    
    return A, rows, cols, essential
    

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