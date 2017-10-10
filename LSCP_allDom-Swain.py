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
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist
from ortools.linear_solver import pywraplp

def RunLSCPCppStyleAPI(optimization_problem_type, SD):
    
    """ Example of simple MCLP program with the C++ style API."""
    solver = pywraplp.Solver('RunIntegerExampleCppStyleAPI', optimization_problem_type)
    
    # Create a global version of:
    # Facility Site Variable X
    X = [None] * numSites
    
    #print sites
    #print np.shape(sites)
    start_time = time.time()
    
    computeCoverageMatrix(SD)
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
    global facilityIDs
    
    # for now, all demands are also sites
    allFD3 = True
    
    # Pull out just the site/demand IDs from the data
    facilityIDs = sites[:,0]
    
    # Pull out just the coordinates from the data
    xyPointArray = sites[:,[1,2]]
    #A = [xyPointArray[i][:] for i in demandIDs]
    #B = [xyPointArray[j][:] for j in facilityIDs]
    A = xyPointArray
    B = A
    #print A
    
    # Compute the distance matrix, using the squared distance
    sqDistMatrix = cdist(A, B,'sqeuclidean')
    print 'Max Point-to-Point Distance = %f' % np.sqrt(np.amax(sqDistMatrix))
    print 'Mean Point-to-Point Distance = %f' % np.sqrt(np.mean(sqDistMatrix))
    print np.shape(sqDistMatrix)
    
    distances = np.unique(sqDistMatrix)
    print np.size(distances)
    
    colmax = np.amax(sqDistMatrix,0)
    minmax = np.amin(colmax)
    
    # print colmax
    print minmax
    
    print "The element in the distances set of the minmax is"
    print np.where(distances==minmax)
    
    print "The site of the minmax is"
    print np.where(colmax==minmax)[0]+1
    
    SDsquared = SD*SD
    TwoSDsquared = 4*SDsquared

    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)

    # Determine neighborhood of sites within SD of sites
    if allFD3 == True:
        SDist = C
    else:
        SDist = (cdist(B, B,'sqeuclidean') <= SDsquared).astype(int)
        

    start_time = time.time()
    C, col_keeps, row_keeps = dominationTrim(C, SDist)
    print 'Domination time = %f' % (time.time()-start_time)

    # shorten the facility data sets
    cols = np.nonzero(col_keeps)[0]
    rows = np.nonzero(row_keeps)[0]
    facilityIDs = [facilityIDs[j] for j in cols]
    numSites = len(facilityIDs)
    print numSites
    numDemands = len(rows)
    print numDemands

    # Convert coverage to sparse matrix
    Nrows,Ncols = np.nonzero(C.astype(bool))
    Nsize = len(Nrows)
#
#     return [SD]
    return 0


def dominationTrim(A, SDist):
    
    rows,cols = A.shape
    r_indices = np.array(range(rows))
    c_indices = np.array(range(cols))
    c_keeps = np.ones(cols)
    r_keeps = np.ones(rows)
    
    # lower triangle of coverage matrix for checking only columns within SD
    # Explanation:
    # looking down each column, each row with a 1 represents a column withing SD of that column
    # using tril means you don't check backwards
    # NOTE: THIS WORKS FOR ONLY SQUARE COVERAGE MATRICES WHERE ALL DEMANDS ARE SITES
    # UPDATE WITH SITE vs. SITE COVERAGE MATRIX FOR OTHER CASES
    L = np.tril(SDist,-1)   # lower triangle matrix
    U = np.triu(SDist,1)    # upper triangle matrix
    
    # start_time = time.time()
    # create a list of sets containing the indices of non-zero elements of each column
    C = csc_matrix(A)
    D = [set(C.indices[C.indptr[i]:C.indptr[i+1]]) for i in range(len(C.indptr)-1)]
    # print 'Matrix to List of Sets CSC Time = %f' % (time.time()-start_time)
    
    # Column domination
    # find subsets, ignoring columns that are known to already be subsets
    for i in range(cols):
        if c_keeps[i]==0:
            continue
        col1 = D[i]
        for j in np.nonzero(L[:,i])[0]:
            # I tried `if keeps[j]==false: continue` here, but that was slower
            # if keeps[j]==False: continue
            col2 = D[j]
            if col2.issubset(col1):
                c_keeps[j] = 0
            elif col1.issubset(col2):
                c_keeps[i] = 0
                break
    E = A[:,c_keeps.astype(bool)]   # delete columns from coverage matrix
    #L = L[:,c_keeps.astype(bool)]   # delete colums from lower triangle matrix
    U = U[:,c_keeps.astype(bool)]   # delete colums from upper triangle matrix
    r_indices = r_indices[c_keeps.astype(bool)]
    
    # Row Domination
    # find subsets, ignoring rows that are known to already be subsets
    # create a list of sets containing the indices of non-zero elements of each column
    R = csr_matrix(E)
    S = [set(R.indices[R.indptr[i]:R.indptr[i+1]]) for i in range(len(R.indptr)-1)]
    
    for i in range(rows):
        if r_keeps[i]==0:
            continue
        row1 = S[i]
        for j in np.nonzero(U[i,:])[0]:
            # I tried `if keeps[j]==false: continue` here, but that was slower
            # if keeps[j]==False: continue
            col2 = S[j]
            if col2.issubset(col1):
                r_keeps[i] = 0
            elif col1.issubset(col2):
                r_keeps[j] = 0
                break
                
    T = E[r_keeps.astype(bool),:]   # delete rows from coverage matrix
    
    rowSumIsOne = np.where(np.sum(T, axis=1)==1)[0]
    nonZeroCols = np.where(T[rowSumIsOne,:])[1]
    print rowSumIsOne
    print nonZeroCols
    print np.unique(nonZeroCols)
    
    return T, c_keeps, r_keeps
    

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
        name = "X,%d" % facilityIDs[j]
        X[j] = solver.BoolVar(name)
        # add the site location variables to the objective function
        objective.SetCoefficient(X[j],1)
    
    # if facility is fixed into the solution, add a constraint to make it so
    for k in range(numForced):
          c3[k] = solver.Constraint(1,1)
          c3[k].SetCoefficient(X[forcedFacilities[k]],1)
    
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
    count = -1
    for j in range(numSites):
        if (X[j].SolutionValue() == 1.0):
            print "Site selected %d" % int(facilityIDs[j])

def read_problem(file):
    global numSites
    global numDemands
    global sites
    global numForced
    
    print 'readFile({0})'.format(file)
    
    lineCount = 0
    i = 0
    numForced = 0
    
    # Use With Statement to automatically close the 'read file' when finished.
    with open(file,'r') as f:
        for line in f:
            line = line.strip()
            
            # ignore comments
            if (line[0] == '#' or line[0] == '%' or len(line) == 0):
                continue
            #print line
            
            if (lineCount == 0):
                # Set the number of sites from the file
                numSites = int(line)
                numDemands = numSites
                
                # Create and instantiate the array 'sites'
                #sites = [[None for k in range(4)] for j in range(numSites)]
                sites = np.empty([numSites,4])
            else:
                row = line.split(" ")
                # Set constraint coefficients
                for j in range(0,4):
                    sites[i,j] = float(row[j])
                i += 1
            lineCount += 1
        # NOTE: CODE BREAKS IF THERE ARE EMPTY LINES AFTER DATA, THIS SHOULD BE FIXED
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
    file = sys.argv[2]
    SD = float(sys.argv[1])
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  elif len(sys.argv) > 1 and len(sys.argv) <= 2:
    SD = float(sys.argv[1])
    file = r'./data/swain.txt'
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  else:
    print "Please Pass: Service Distance; Data to Use"
    print "Problem not executed!"