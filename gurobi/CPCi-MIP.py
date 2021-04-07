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

# this sets the number of concurrency and threads for parallel computation
# default is 0
threads = 0
conc = 0

#threads = 1
#setParam(GRB.Param.Threads, threads)
#conc = 0
#setParam(GRB.Param.ConcurrentMIP, conc)


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
    solution[p-1,1] = SDmin
    
    BuildModel(m, 0, distMatrix)
    print('  p, SD')
    displaySolution(p, solution[p-1,1])
    
    for p in range(2, numSites):

        # update the right hand side of the facility constraint
        m.getConstrByName("c1").setAttr(GRB.Attr.RHS, p)
        
        SolveModel(m)
        SDmin = m.objVal
        solution[p-1,1] = SDmin
        
        displaySolution(p, solution[p-1,1])
    
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
    #A = [xyPointArray[i][:] for i in demandIDs]
    #B = [xyPointArray[j][:] for j in siteIDs]
    A = xyPointArray
    B = A
    #print A
    
    # Compute the distance matrix, using the euclidean distance
    distMatrix = np.ceil(cdist(A, B,'euclidean')).astype(int)
    return distMatrix

def BuildModel(m, p, d):
    
    rNumD = range(numDemands)
    rNumS = range(numSites)
    
    # DECLARE VARIABLES:
    # Assignment variables X
    # =1 if demand i is assigned to facility j
    X = m.addVars(numDemands, numSites,
                  vtype=GRB.BINARY)
                  
    # Facility Site binary decision variables Y
    # =1 if facility is located at site j
    Y = m.addVars(numSites,
                  vtype=GRB.BINARY)

    # Cover distance variable Z
    # continuous variable to be minimized
    Z = m.addVar(vtype=GRB.CONTINUOUS, obj = 1.0)
    
    # Define Facility Constraint (c1):
    #m.addConstr(quicksum(Y[j] for j in range(numSites)) <= p, "c1") # uses old notation style
    m.addConstr(Y.sum() <= p, "c1")   # uses new tupledict notation style
    
    ### FYI: Making constraints c2 and c4 using Tupledict notation makes the model slower
    # Define Assignment Constraints (c2)
    # Define Z to be the largest distance from any demand to any facility (c4)
    
    # make constraint 2 using tupledict notation
    # m.addConstrs((X.sum(i, '*') == 1 for i in rNumD), "c2")
    #
    # # make constraint 4 using tupledict notation
    # dDict = dict(((i,j), d[i][j]) for i, j in itertools.product(rNumD, rNumS))
    # m.addConstrs((X.prod(dDict, i, '*') <= Z for i in rNumD), "c4")

    for i in rNumD:
        m.addConstr(quicksum(X[i,j] for j in rNumS) == 1, "c2[%d]" % i)
        m.addConstr(quicksum(X[i,j]*d[i,j] for j in rNumS) - Z <= 0, "c4[%d]" % i)

        for j in rNumS:
            # add the balinsky assignment constraints (c3)
            # Yj - Xij >= 0 <--- canonical form of the assignment constraint
            m.addConstr(X[i,j] <= Y[j], "c3[%d,%d]" % (i,j))

    # The objective is to minimize the number of located facilities
    m.modelSense = GRB.MINIMIZE

    m.update()
    print('Number of variables = %d' % m.numvars)
    print('Number of constraints = %d' % m.numconstrs)
    #m.printStats()
    
    print
    return 0

def SolveModel(m):
    """Solve the problem and print the solution."""
    # m.Params.ResultFile = "output.sol"
    m.optimize()
    
    
def displaySolution(p, SDmin):
    # The objective value and the minimum service distance
    print('%3d, %d' % (p, SDmin))
    

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
    
    print('%d locations, %d threads'% (numSites,threads))
    print()


def main(unused_argv):
    print('---- CPC-MIP with Gurobi -----')
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