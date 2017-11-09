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

def RunLSCPCppStyleAPI(optimization_problem_type):
    
    """ Example of simple MCLP program with the C++ style API."""
    solver = pywraplp.Solver('RunIntegerExampleCppStyleAPI', optimization_problem_type)
    
    start_time = time.time()
    
    sqDistances, sqDistMatrix = computeDistances()
    
    print '  p, SD'   
    print '%3d, 0' % numSites
    
    solution = np.empty([numSites, 2])
    solution[:,0] = range(numSites, 0, -1)
    solution[0,1] = 0
    currP = numSites
    iters = 0
    #print solution
    
    for i in range(1,len(sqDistances)):
        SDsquared = sqDistances[i]
        computeCoverageMatrix(sqDistMatrix, SDsquared)
    
        # Facility Site Variable X
        X = [None] * numSites

        BuildModel(solver, X)
        SolveModel(solver)

        # get the solution and clear the solver
        p = solver.Objective().Value()
        solver.Clear()
        
        # check the output
        while (p < currP):
            currP -= 1
            displaySolution(currP, SDsquared)

        # terminate the search when p == 1
        if (p == 2):
            p = 1
            SDsquared = np.amin(np.amax(sqDistMatrix,0))
            displaySolution(p, SDsquared)
            iters = i+1
            break
        
    total_time = time.time()-start_time
    print
    print '%d LSCP distances evaluated' % iters
    print 'Total problem solved in %f seconds' % total_time
    print
    
def computeDistances():
        
    #declare a couple variables
    global distances
    global Nrows
    global Ncols
    global Nsize
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
    print sqDistMatrix.shape
    
    # print 'Max Point-to-Point Distance = %f' % np.sqrt(np.amax(sqDistMatrix))
    # print 'Mean Point-to-Point Distance = %f' % np.sqrt(np.mean(sqDistMatrix))
    # print np.shape(sqDistMatrix)
    #
    sqDistances = np.unique(sqDistMatrix)
    
    # print np.size(distances)
    #
    # colmax = np.amax(sqDistMatrix,0)
    # minmax = np.amin(colmax)
    #
    # # print colmax
    # print minmax**(0.5)
    #
    # print "The element in the distances set of the minmax is"
    # print np.where(distances==minmax)
    #
    # print "The site of the minmax is"
    # print np.where(colmax==minmax)[0]+1
    
    return sqDistances, sqDistMatrix
    

def computeCoverageMatrix(sqDistMatrix, SDsquared):
        
    #declare a couple variables
    global Nrows
    global Ncols
    global Nsize
    
    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)
        
    # Convert coverage to sparse matrix
    Nrows,Ncols = np.nonzero(C.astype(bool))
    Nsize = len(Nrows)

    return 0

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
    
    # print 'Number of variables = %d' % solver.NumVariables()
    # print 'Number of constraints = %d' % solver.NumConstraints()
    # print
    return 0

def SolveModel(solver):
    """Solve the problem and print the solution."""
    result_status = solver.Solve()

    # The problem has an optimal solution.
    assert result_status == pywraplp.Solver.OPTIMAL, "The problem does not have an optimal solution!"

    # The solution looks legit (when using solvers others than
    # GLOP_LINEAR_PROGRAMMING, verifying the solution is highly recommended!).
    assert solver.VerifySolution(1e-7, True)
    
def displaySolution(p, SDsquared):

    # The objective value and the minimum service distance
    print '%3d, %f' % (p, SDsquared**0.5)
    # print the selected sites
    

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
    

def Announce(solver, api_type):
    print ('---- P-Center LSCP with ' + solver + ' (' +
        api_type + ') -----')

def RunSCIP_LSCPexampleCppStyleAPI():
    if hasattr(pywraplp.Solver, 'SCIP_MIXED_INTEGER_PROGRAMMING'):
        Announce('SCIP', 'C++ style API')
        RunLSCPCppStyleAPI(pywraplp.Solver.SCIP_MIXED_INTEGER_PROGRAMMING)

def RunCBC_LSCPexampleCppStyleAPI():
    if hasattr(pywraplp.Solver, 'CBC_MIXED_INTEGER_PROGRAMMING'):
        Announce('CBC', 'C++ style API')
        RunLSCPCppStyleAPI(pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

def RunBOP_LSCPexampleCppStyleAPI():
    if hasattr(pywraplp.Solver, 'BOP_INTEGER_PROGRAMMING'):
        Announce('BOP', 'C++ style API')
        RunLSCPCppStyleAPI(pywraplp.Solver.BOP_INTEGER_PROGRAMMING)


def main(unused_argv):
    RunCBC_LSCPexampleCppStyleAPI()
    #RunSCIP_LSCPexampleCppStyleAPI()
    #RunBOP_LSCPexampleCppStyleAPI()


""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 1 and len(sys.argv) <= 2:
    file = './data/' + sys.argv[1]
    print
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  elif len(sys.argv) > 0 and len(sys.argv) <= 1:
    file = './data/swain.dat'
    print
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  else:
    print "Please Pass: Service Distance; Data to Use"
    print "Problem not executed!"