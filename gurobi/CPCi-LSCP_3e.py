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
import pCenterBrute as brute
from scipy.spatial.distance import cdist
from gurobipy import *
setParam('OutputFlag', 0)   # mute solver meta-info

def Run_pCenterLSCP():
    
    """Example of complete p-Center program with the Gurobi API"""
    m = Model()
    
    start_time = time.time()
    
    distances, distMatrix = computeDistances()
    
    # p = numSites, SD = 0 is a trivial solution
    print('  p, SD')
    p = numSites
    SD = 0   
    displaySolution(p, SD)
    
    solution = np.empty([numSites, 2])
    # solution[:,0] = range(1, numSites+1)
    # solution[p-1,1] = 0
    currP = numSites
    
    SD = distances[0]
    C = computeCoverageMatrix(distMatrix, SD)

    BuildModel(m)
    SolveModel(m)
    p = m.objVal
    
    while (p < currP):
        currP -= 1
        solution[currP-1,1] = SD
        displaySolution(currP, SD)

    for k in range(1,len(distances)):
        SD = distances[k]

        diff, C = updateCoverCoefficeints(distMatrix, SD, C)
                
        for i in range(numDemands):
            for j in diff[i]:
                m.chgCoeff(m.getConstrByName("c[%d]" % i), X[j], 1)
        
        SolveModel(m)

        # get the solution and clear the solver
        p = m.objVal

        # check the output
        while (p < currP):
            currP -= 1
            solution[currP-1,1] = SD
            displaySolution(currP, SD)

        # solve brute force for p == 3
        if (p == 4):
            p = 3
            if numSites > 263:
                SD, rows = brute.nbParallel3(distMatrix, numSites)
            elif numSites > 132:
                SD, rows = brute.nbSerial3(distMatrix, numSites)
            else:
                SD, rows = brute.chunk3(distMatrix, numSites)
            
            displaySolution(p, SD)
            iters = k+1

        # solve brute force for p == 2
        if (p == 3):
            p = 2
            if numSites > 1100:
                SD, rows = brute.nbParallel2(distMatrix, numSites)
            elif numSites > 507:
                SD, rows = brute.nbSerial2(distMatrix, numSites)
            else:
                SD, rows = brute.chunk2(distMatrix, numSites)
            
            displaySolution(p, SD)
            iters = k+1

        # solve brute force for p == 1
        if (p == 2):
            p = 1
            SD = np.amin(np.amax(distMatrix,0))
            solution[p-1,1] = SD
            displaySolution(p, SD)
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
    
def computeDistances():
    
    # Pull out just the site/demand IDs from the data
    siteIDs = sites[:,0]
    
    # Pull out just the coordinates from the data
    xyPointArray = sites[:,[1,2]]

    A = xyPointArray
    B = A
    
    # Compute the distance matrix, using the squared distance
    distMatrix = np.ceil(cdist(A, B,'euclidean')).astype(int)
    sqDistances = np.unique(distMatrix)
    
    return sqDistances, distMatrix
    

def computeCoverageMatrix(distMatrix, SD):
        
    #declare a couple variables
    global cover_rows
    
    # Determine neighborhood of demands within SD of sites
    C = (distMatrix <= SD).astype(int)
        
    # Convert coverage to array of nonzero elements in each row
    cover_rows = [np.nonzero(t)[0] for t in C]

    return C


def updateCoverCoefficeints(distMatrix, SD, B):
    
    # Determine neighborhood of demands within SD of sites
    C = (distMatrix <= SD).astype(int)
    diff = [np.nonzero(t)[0] for t in (C-B)]
    
    return diff, C


def BuildModel(m):
    
    global X
    
    # DECLARE VARIABLES:
    # Facility Site binary decision variables X
    # Each has a coefficient of 1 in the objective
    X = m.addVars(numSites,
                  vtype=GRB.BINARY,
                  obj=1)
    
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
    
def displaySolution(p, SD):
    print('%3d. %d' % (p, SD))

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


""" Main will take in 2 argument: Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 1 and len(sys.argv) <= 2:
    file = '../data/' + sys.argv[1]
    print()
    print("Problem instance from: ", file)
    read_problem(file)
    main(sys.argv[1])
  elif len(sys.argv) > 0 and len(sys.argv) <= 1:
    file = '../data/swain.dat'
    print()
    print("Problem instance from: ", file)
    read_problem(file)
    main('swain.dat')
  else:
    print("Please Pass: Data to Use")
    print("Problem not executed!")