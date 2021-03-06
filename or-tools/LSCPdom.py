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
from scipy.spatial.distance import cdist
from ortools.linear_solver import pywraplp

def RunLSCPCppStyleAPI(optimization_problem_type, SD):
    
    """ Example of simple MCLP program with the C++ style API."""
    solver = pywraplp.Solver('RunIntegerExampleCppStyleAPI', optimization_problem_type)
    
    #print sites
    #print np.shape(sites)
    start_time = time.time()
    
    computeCoverageMatrix(SD)

    # Facility Site Variable X
    X = [None] * numSites
    
    BuildModel(solver, X)
    SolveModel(solver)
    
    total_time = time.time()-start_time
    p = solver.Objective().Value()
    
    displaySolution(X, p, total_time)
    
    
    
def computeCoverageMatrix(SD):
        
    #declare a couple variables
    global distances
    global numDemands
    global numSites
    global Nrows
    global Ncols
    global Nsize
    global cols
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
    sqDistMatrix = cdist(A, B,'sqeuclidean')
    distances = np.unique(sqDistMatrix)
    SDsquared = SD*SD

    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)
    # Determine neighborhood of sites within 2*SD of sites (symmetric)
    C2 = (sqDistMatrix <= 4*SDsquared).astype(int)
        

    start_time = time.time()
    C, cols = dominationTrim(C, C2)
    print 'Domination time = %f' % (time.time()-start_time)

    # shorten the facility data sets
    siteIDs = [siteIDs[j] for j in cols]
    numSites = len(siteIDs)

    # Convert coverage to sparse matrix
    Nrows,Ncols = np.nonzero(C.astype(bool))
    Nsize = len(Nrows)

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
    

def BuildModel(solver, X):
    
    infinity = solver.infinity()
    
    # DECLARE CONSTRAINTS:
    # declare demand coverage constraints (binary integer: 1 if UNCOVERED, 0 if COVERED)
    c1 = [None]*numDemands
    
    # declare the objective
    objective = solver.Objective()
    objective.SetMinimization()
    
    # initialize the X variables as Binary Integer (Boolean) variables
    for j in range(numSites):
        name = "X,%d" % siteIDs[j]
        X[j] = solver.BoolVar(name)
        # add the site location variables to the objective function
        objective.SetCoefficient(X[j],1)
    
    # add demands to the objective and coverage constraints
    for i in range(numDemands):
        # Covering constraints
        c1[i] = solver.Constraint(1, solver.infinity())

    # add facility coverages to the coverage constraints
    for k in range(Nsize):
        c1[Nrows[k]].SetCoefficient(X[Ncols[k]],1)
    
    print 'Number of variables = %d' % solver.NumVariables()
    print 'Number of constraints = %d' % solver.NumConstraints()
    print
    return 0

def SolveModel(solver):
    """Solve the problem and print the solution."""
    result_status = solver.Solve()

    # The problem has an optimal solution.
    assert result_status == pywraplp.Solver.OPTIMAL, "The problem does not have an optimal solution!"

    # The solution looks legit (when using solvers others than
    # GLOP_LINEAR_PROGRAMMING, verifying the solution is highly recommended!).
    assert solver.VerifySolution(1e-7, True)
    
def displaySolution(X, p, total_time):

    print 'Total problem solved in %f seconds' % total_time
    print
    # The objective value of the solution.
    print 'p = %d' % p
    print 'SD = %f' % SD
    # print the selected sites
    print
    for j in range(numSites):
        if (X[j].SolutionValue() == 1.0):
            print "Site selected %d" % int(siteIDs[j])
            
    # plot solution
    plot.plotSolution(sites, X, cols, SD)
    

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


def Announce(solver, api_type):
    print ('---- Integer programming example with ' + solver + ' (' +
        api_type + ') -----')

def RunSCIP_LSCPExampleCppStyleAPI(SD):
    if hasattr(pywraplp.Solver, 'SCIP_MIXED_INTEGER_PROGRAMMING'):
        Announce('SCIP', 'C++ style API')
        RunLSCPCppStyleAPI(pywraplp.Solver.SCIP_MIXED_INTEGER_PROGRAMMING, SD)

def RunCBC_LSCPexampleCppStyleAPI(SD):
    if hasattr(pywraplp.Solver, 'CBC_MIXED_INTEGER_PROGRAMMING'):
        Announce('CBC', 'C++ style API')
        RunLSCPCppStyleAPI(pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING, SD)

def RunBOP_LSCPexampleCppStyleAPI(SD):
    if hasattr(pywraplp.Solver, 'BOP_INTEGER_PROGRAMMING'):
        Announce('BOP', 'C++ style API')
        RunLSCPCppStyleAPI(pywraplp.Solver.BOP_INTEGER_PROGRAMMING, SD)


def main(unused_argv):
    RunCBC_LSCPexampleCppStyleAPI(SD)
    #RunSCIP_LSCPexampleCppStyleAPI(SD)
    #RunBOP_LSCPexampleCppStyleAPI(SD)


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