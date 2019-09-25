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
import readDataFiles
import plot
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
    displaySolution(p, SDsquared, 0)
    
    solution = np.empty([numSites, 2])
    start_time_mini = time.time()
    # solution[:,0] = range(1, numSites+1)
    # solution[p-1,1] = 0
    currP = numSites
    
    SDsquared = sqDistances[0]
    C = computeCoverageMatrix(sqDistMatrix, SDsquared)

    BuildModel(m)
    SolveModel(m)
    p = m.objVal
    
    while (p < currP):
        currP -= 1
        total_time_mini = time.time()-start_time_mini
        solution[currP-1,1] = SDsquared**0.5
        displaySolution(currP, SDsquared, total_time_mini)
        start_time_mini = time.time()

    for k in range(1,len(sqDistances)):
        SDsquared = sqDistances[k]

        diff, C = updateCoverCoefficeints(sqDistMatrix, SDsquared, C)
                
        for i in range(numDemands):
            for j in diff[i]:
                m.chgCoeff(m.getConstrByName("c[%d]" % i), X[j], 1)
        
        SolveModel(m)

        # get the solution and clear the solver
        p = m.objVal

        # check the output
        while (p < currP):
            currP -= 1
            total_time_mini = time.time()-start_time_mini
            solution[currP-1,1] = SDsquared**0.5
            displaySolution(currP, SDsquared, total_time_mini)
            start_time_mini = time.time()
                
        # terminate the search when p == 1
        if (p == 2):
            p = 1
            total_time_mini = time.time()-start_time_mini
            SDsquared = np.amin(np.amax(sqDistMatrix,0))
            solution[p-1,1] = SDsquared**0.5
            displaySolution(p, SDsquared, total_time_mini)
            iters = k+1
            break
        if (p == 1):
            iters = k
            break
        
    total_time = time.time()-start_time
    #print solution
    print
    print('%d LSCP distances evaluated' % iters)
    print('Total problem solved in %f seconds' % total_time)
    print()
    #plot.plotTradeoff(file, solution)
    
def computeDistances():
        
    #declare a couple variables
    global distances
    global siteIDs
    
    # Pull out just the site/demand IDs from the data
    siteIDs = sites[:,0]
    
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
    global cover_rows
    
    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)
        
    # Convert coverage to array of nonzero elements in each row
    cover_rows = [np.nonzero(t)[0] for t in C]

    return C


def updateCoverCoefficeints(sqDistMatrix, SDsquared, B):
    
    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)
    diff = [np.nonzero(t)[0] for t in (C-B)]
    
    return diff, C


def BuildModel(m):
    
    global X
    
    # DECLARE VARIABLES:
    # Facility Site binary decision variables X
    # Each has a coefficient of 1 in the objective
    X = m.addVars(numSites,
                  vtype=GRB.BINARY,
                  obj=np.ones(numSites),
                  name="X")
    
    # Define Coverage Constraints:
    for i in range(numDemands):
        m.addConstr(quicksum(X[j] for j in cover_rows[i]) >= 1, "c[%d]" % i)
    
    # The objective is to minimize the number of located facilities
    m.modelSense = GRB.MINIMIZE
    
    # m.update()
    # print 'Number of variables = %d' % solver.NumVariables()
    # print 'Number of constraints = %d' % solver.NumConstraints()
    # print
    return 0


def SolveModel(m):
    """Solve the problem and print the solution."""
    # m.Params.ResultFile = "output.sol"
    m.optimize()
    
    
def displaySolution(p, SDsquared, time):
    # The objective value and the minimum service distance
    print('%3d, %f, %f' % (p, SDsquared**0.5, time))
    

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
    print('---- pCenterLSCP with Gurobi -----')
    Run_pCenterLSCP()


""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 1 and len(sys.argv) <= 2:
    file = '../data/' + sys.argv[1]
    print
    print("Problem instance from: ", file)
    read_problem(file)
    main(sys.argv[1])
  elif len(sys.argv) > 0 and len(sys.argv) <= 1:
    file = '../data/swain.dat'
    print
    print("Problem instance from: ", file)
    read_problem(file)
    main('swain.dat')
  else:
    print("Please Pass: Service Distance; Data to Use")
    print("Problem not executed!")