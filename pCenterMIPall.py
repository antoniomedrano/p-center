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

def RunMIPCppStyleAPI(optimization_problem_type, p):
    
    """ Example of simple MCLP program with the C++ style API."""
    solver = pywraplp.Solver('RunIntegerExampleCppStyleAPI', optimization_problem_type)
    
    start_time = time.time()
    
    distMatrix = computeCoverageMatrix()

    # Facility Site Variable X
    X = [[None for j in range(numSites)] for i in range(numSites)]
    Y = [None] * numSites
    Z = None

    BuildModel(solver, X, Y, Z, p, distMatrix)
    SolveModel(solver)
    
    total_time = time.time()-start_time
    SDmin = solver.Objective().Value()
    
    print
    print 'Total problem solved in %f seconds' % total_time
    print
    #displaySolution(Y, p, SDmin, total_time)
    
    
def computeCoverageMatrix():
        
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
    distMatrix = cdist(A, B,'euclidean')
    #
    # distances = np.unique(distMatrix)
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
    
    
    # # Determine neighborhood of demands within SD of sites
    # C = (sqDistMatrix <= SDsquared).astype(int)
    #
    # # Convert coverage to sparse matrix
    # Nrows,Ncols = np.nonzero(C.astype(bool))
    # Nsize = len(Nrows)

    return distMatrix

def BuildModel(solver, X, Y, Z, p, d):
    
    infinity = solver.infinity()
    
    # DECLARE CONSTRAINTS:
    # declare demand coverage constraints (binary integer: 1 if UNCOVERED, 0 if COVERED)
    c1 = None  # number of facilities constraint
    c2 = [None] * numSites  # assign sites to one facility
    c3 = [None] * numSites**2  # assignment only to located facilities
    c4 = [None] * numSites  # force Z to be > the distance from any client to the assigned facility
    
    # declare the objective
    objective = solver.Objective()
    objective.SetMinimization()
    # declare the Z linear variable and assign it to the objective
    Z = solver.NumVar(0, infinity, 'Z')
    objective.SetCoefficient(Z, 1)
    
    # <= constraint for locating p facilities
    c1 = solver.Constraint(0,p)
    
    for j in range(numSites):
        # initialize the Y facility location variables
        name = "Y,%d" % siteIDs[j]
        Y[j] = solver.BoolVar(name)
        # set coefficients for Y variables in constraint 1
        c1.SetCoefficient(Y[j],1)
        
    # initialize the X variables as Binary Integer (Boolean) variables
    for i in range(numSites):
        # Covering constraints = 1
        c2[i] = solver.Constraint(1, 1)
        # Z distance assignment constraint >= 0
        c4[i] = solver.Constraint(0, infinity)
        c4[i].SetCoefficient(Z, 1)
        
        for j in range(numSites):
            # initialize the Xij assignment variables
            name = "X,%d,%d" % (siteIDs[i], siteIDs[j])
            X[i][j] = solver.BoolVar(name)

            # set the variable coefficients of the sum(Xij) = 1 for each i
            c2[i].SetCoefficient(X[i][j], 1)
           
            # add the balinsky assignment constraints
            # Yj - Xij >= 0 <--- canonical form of the assignment constraint
            c3[i*numSites+j] = solver.Constraint(0, infinity) # c2 rhs
            c3[i*numSites+j].SetCoefficient(X[i][j], -1)
            c3[i*numSites+j].SetCoefficient(Y[j], 1)
            
            # add the assignment distance variable coefficients
            c4[i].SetCoefficient(X[i][j], -d[i,j])
    
    # # add facility coverages to the coverage constraints
    # for k in range(Nsize):
    #     c1[Nrows[k]].SetCoefficient(X[Ncols[k]],1)
    
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
    
def displaySolution(p, SDmin):

    # The objective value and the minimum service distance
    print '%3d, %f' % (p, SDmin)
    
    print 'Total problem solved in %f seconds' % total_time
    print
    # The objective value of the solution.
    print 'p = %d' % p
    print 'SD = %f' % SDmin
    # print the selected sites
    print    
    for j in range(numSites):
        if (Y[j].SolutionValue() == 1.0):
            print "Site selected %d" % int(siteIDs[j])
    
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
        print 'Error reading file'
        raise
        
    numSites = sites.shape[0]    
    numDemands = numSites
    
    #plot.plotData(sites)
    
    print '%d locations' % numSites
    print 'Finished Reading File!'


def Announce(solver, api_type):
    print ('---- P-Center MIP with ' + solver + ' (' +
        api_type + ') -----')

def RunSCIP_MIPexampleCppStyleAPI():
    if hasattr(pywraplp.Solver, 'SCIP_MIXED_INTEGER_PROGRAMMING'):
        Announce('SCIP', 'C++ style API')
        RunMIPCppStyleAPI(pywraplp.Solver.SCIP_MIXED_INTEGER_PROGRAMMING, )

def RunCBC_MIPexampleCppStyleAPI():
    if hasattr(pywraplp.Solver, 'CBC_MIXED_INTEGER_PROGRAMMING'):
        Announce('CBC', 'C++ style API')
        RunMIPCppStyleAPI(pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING, )

def RunBOP_MIPexampleCppStyleAPI():
    if hasattr(pywraplp.Solver, 'BOP_INTEGER_PROGRAMMING'):
        Announce('BOP', 'C++ style API')
        RunMIPCppStyleAPI(pywraplp.Solver.BOP_INTEGER_PROGRAMMING, )


def main(unused_argv):
    RunCBC_MIPexampleCppStyleAPI()
    #RunSCIP_MIPexampleCppStyleAPI()
    #RunBOP_MIPexampleCppStyleAPI()


""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 1 and len(sys.argv) <= 2:
    file = './data/' + sys.argv[1]
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  elif len(sys.argv) > 0 and len(sys.argv) <= 1:
    p = float(sys.argv[1])
    file = './data/swain.dat'
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  else:
    print "Please Pass: Service Distance; Data to Use"
    print "Problem not executed!"