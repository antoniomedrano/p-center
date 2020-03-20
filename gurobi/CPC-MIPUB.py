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
from scipy.spatial.distance import cdist
from gurobipy import *
setParam('OutputFlag', 0)   # mute solver meta-info

threads = 0
conc = 0

# threads = 6
# setParam(GRB.Param.Threads, threads)
# conc = 0
# setParam(GRB.Param.ConcurrentMIP, conc)   

def Run_pCenter():
    
    """Example of complete p-Center program with the Gurobi API"""
    m = Model()
    
    start_time = time.time()
    
    distMatrix = computeDistanceMatrix()

    solution = np.empty([numSites, 2])
    solution[:,0] = range(1, numSites+1)

    # p = 1 is a trivial solution min(max(dist))
    p = 1
    SDmin = np.amin(np.amax(distMatrix,0))
    solution[p-1,1] = SDmin**0.5

    C = computeCoverageMatrix(distMatrix, SDmin)
    BuildModel(m, 2, distMatrix)
    
    print('  p, SD')
    displaySolution(p, solution[p-1,1])
    
    p = 2
    SolveModel(m)
    SDmin = m.objVal
    solution[p-1,1] = SDmin**0.5
    displaySolution(p, solution[p-1,1])
    
    for i in range(3, numSites):
        p = i
        
        # find the difference in the coverage matrix from p=i+1 to p=i
        # add a small amount to avoid issues with numerical truncation
        diff, C = updateCoverCoefficeints(distMatrix, SDmin+.00001, C)

        # for u574 and pr439, I had to add 0.0005 to SDmin so it wouldn't crash.
        # 0.0001 crashes at p=282, could maybe add a bit more just then but I don't like that.
        # diff, C = updateCoverCoefficeints(distMatrix, SDmin + 0.0002, C) #u574 maybe add less
        # diff, C = updateCoverCoefficeints(distMatrix, SDmin + 0.02, C) #pr439
 
        # query the variables from the previous solution to perform an explicit warm start
        for v in m.getVars():
            v._prev = v.X
        
        # update the right hand side of the facility constraint
        m.getConstrByName("c1").setAttr(GRB.Attr.RHS, p)
        for i in range(numDemands):
            # for assignment variables associated with d > SDmin
            for j in diff[i]:
                # delete the variable from the model
                m.remove(m.getVarByName("x[%d,%d]" % (i,j)))

                # remove it's associated balinski constraint
                m.remove(m.getConstrByName("c3[%d,%d]" % (i,j)))

        m.update()
        # assign the starting variable values from the previous solution to do a warm start
        for v in m.getVars():
            v.Start = v._prev
        
        SolveModel(m)
        SDmin = m.objVal
        solution[p-1,1] = SDmin**0.5
        
        displaySolution(p, SDmin**0.5)
    
    # solution for p = numSites is SDmin = 0    
    solution[numSites-1,1] = 0
    displaySolution(numSites, 0)
        
    total_time = time.time()-start_time
    
    print
    print('Problem solved in %f secs with %d threads and concurrency of %d' % (total_time, threads, conc))
    print
    
    
def computeDistanceMatrix():
    
    # Pull out just the coordinates from the data
    xyPointArray = sites[:,[1,2]]
    A = xyPointArray
    B = A
    #print A
    
    # Compute the distance matrix, using the euclidean distance
    distMatrix = cdist(A, B,'sqeuclidean')
    
    return distMatrix
    
    
def computeCoverageMatrix(distMatrix, SD_UB):
    
    global cover_rows

    # Determine neighborhood of demands within SD of sites
    C = (distMatrix <= SD_UB).astype(int)

    # Convert coverage to sparse matrix
    cover_rows = [np.nonzero(t)[0] for t in C]
    
    return C
    
def updateCoverCoefficeints(distMatrix, SD_UB, B):
    
    # Determine neighborhood of demands within SD of sites
    C = (distMatrix <= SD_UB).astype(int)
    diff = [np.nonzero(t)[0] for t in (C-B)]
    
    return diff, C


def BuildModel(m, p, d):
        
    # DECLARE VARIABLES:
    # Assignment variables X
    # =1 if demand i is assigned to facility j
    X = {}
    for i in range(numDemands):
        for j in cover_rows[i]:
            X[i,j] = m.addVar(vtype=GRB.BINARY, name="x[%d,%d]" % (i,j))
    
    # Facility Site binary decision variables Y
    # =1 if facility is located at site j
    Y = m.addVars(numSites,
                  vtype=GRB.BINARY)

    # Cover distance variable Z
    # continuous variable to be minimized
    Z = m.addVar(vtype=GRB.CONTINUOUS, obj = 1.0)
    
    # Define Facility Constraint (c1):
    # m.addConstr(Y[j] for j in range(numSites)) <= p, "c1") # uses old notation style
    m.addConstr(Y.sum() <= p, "c1")   # uses new tupledict notation style

    # Define Assignment Constraints (c2)
    # Define Z to be the largest distance from any demand to any facility (c4)
    for i in range(numDemands):
        m.addConstr(quicksum(X[i,j] for j in cover_rows[i]) == 1, "c2[%d]" % i)
        m.addConstr(quicksum(X[i,j]*d[i,j] for j in cover_rows[i]) - Z <= 0, "c4[%d]" % i)

        #for j in range(numSites):
        for j in cover_rows[i]:
            # add the balinsky assignment constraints (c3)
            # Yj - Xij >= 0 <--- canonical form of the assignment constraint
            m.addConstr(X[i,j] <= Y[j], "c3[%d,%d]" % (i,j))

    # The objective is to minimize the number of located facilities
    m.modelSense = GRB.MINIMIZE
    
    m.update()
    print('Number of variables = %d' % m.numvars)
    print('Number of constraints = %d' % m.numconstrs)
    #m.printStats()
    
    print()
    return 0

def SolveModel(m):
    """Solve the problem and print the solution."""
    # m.Params.ResultFile = "output.sol"
    m.optimize()
    
    
def displaySolution(p, SDmin):
    # The objective value and the minimum service distance
    print('%3d, %f' % (p, SDmin))
    

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
    print()


def main(unused_argv):
    print('---- CPC-MIPUB with Gurobi -----')
    Run_pCenter()


""" Main will take in 1 argument: Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 1 and len(sys.argv) <= 2:
    file = '../data/' + sys.argv[1]
    print("Problem instance from: ", file)
    read_problem(file)
    main(None)
  elif len(sys.argv) > 0 and len(sys.argv) <= 1:
    file = '../data/swain.dat'
    print("Problem instance from: ", file)
    read_problem(file)
    main(None)
  else:
    print("Please Pass: Data to Use")
    print("Problem not executed!")