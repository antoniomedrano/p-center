# Copyright 2021 Antonio Medrano
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
# Author: F. Antonio Medrano

import sys
import time
import numpy as np
import readDataFiles
import plot
from scipy.spatial.distance import cdist
from gurobipy import *
setParam('OutputFlag', 0)   # mute solver meta-info

# this sets the number of threads for parallel computation
threads = 0
conc = 0
# setParam(GRB.Param.Threads, threads)
# threads = 0
# setParam(GRB.Param.ConcurrentMIP, conc)
# conc = 0
# setParam(GRB.Param.MIPFocus, 1)
# setParam(GRB.Param.MIPGap, 0.01)

def Run_pCenter(p):
    
    """ Example of simple p-Center program with the Gurobi Python API"""
    m = Model()    
    # m.Params.SolutionLimit = 4
    
    start_time = time.time()
    
    distMatrix = computeDistanceMatrix()
    #print distMatrix

    BuildModel(m, p, distMatrix)

    SolveModel(m)
    
    total_time = time.time()-start_time
    
    displaySolution(m, p, total_time)
    
    
def computeDistanceMatrix():
        
    #declare a couple variables
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
    distMatrix = cdist(A, B,'sqeuclidean')

    return distMatrix

def BuildModel(m, p, d):
    
    # DECLARE VARIABLES:
    # Assignment variables X
    # =1 if demand i is assigned to facility j
    sitesRange = range(numSites)
    demandsRange = range(numDemands)
    
    X = m.addVars(numDemands, numSites,
                  vtype=GRB.BINARY,
                  name="X")
    # Facility Site binary decision variables Y
    # =1 if facility is located at site j
    Y = m.addVars(sitesRange,
                  vtype=GRB.BINARY,
                  name="Y")

    # Cover distance variable Z
    # continuous variable to be minimized
    Z = m.addVar(vtype=GRB.CONTINUOUS, obj = 1.0, name="Z")
    
    # Define Facility Constraint (c1):
    m.addConstr(Y.sum() <= p, "c1")   # uses new tupledict notation style

    # Define Assignment Constraints (c2)
    # Define Z to be the largest distance from any demand to any facility (c4)
    for i in range(numDemands):
        m.addConstr(quicksum(X[i,j] for j in range(numSites)) == 1, "c2[%d]" % i)
        m.addConstr(quicksum(X[i,j]*d[i,j] for j in range(numSites)) - Z <= 0, "c4[%d]" % i)

        for j in range(numSites):
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
    
def displaySolution(m, p, total_time):

    print('Problem solved in %f secs with %d threads and concurrency of %d' % (total_time, threads, conc))
    print()
    # The objective value of the solution.
    print('p = %d' % p)
    print('SD = %f' % m.objVal**0.5)
    # print the selected sites
    print()
    for j in range(numSites):
        v = m.getVarByName("Y[%d]" % j)
        if (v.x == 1.0):
            print("Site selected %s" % int(siteIDs[j]))
    
    #plot solution
    # needs fixing for gurobi
    #plot.plotSolution(sites, m.Y, range(numSites), m.objVal**0.5)
    

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
    print('Finished Reading File!')


def main(unused_argv):
    print('---- P-Center with Gurobi -----')
    Run_pCenter(p)


""" Main will take in 2 arguments: p-Facilities; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 2 and len(sys.argv) <= 3:
    file = '../data/' + sys.argv[2]
    p = float(sys.argv[1])
    print("Problem instance from: ", file)
    read_problem(file)
    main(None)
  elif len(sys.argv) > 1 and len(sys.argv) <= 2:
    p = float(sys.argv[1])
    file = '../data/swain.dat'
    print("Problem instance from: ", file)
    read_problem(file)
    main(None)
  else:
    print("Please Pass: number of facilities; Data to Use")
    print("Problem not executed!")