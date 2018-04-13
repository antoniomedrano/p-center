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
from scipy.spatial.distance import cdist
from gurobipy import *

def RunLSCP(SD):
    
    # Example of simple LSCP program with the C++ style API.
    m = Model()
    
    start_time = time.time()
    
    computeCoverageMatrix(SD)

    BuildModel(m)
    SolveModel(m)

    total_time = time.time()-start_time

    p = m.objVal
    displaySolution(m, p, total_time)
    
    
def computeCoverageMatrix(SD):
        
    #declare a couple variables
    global cover_rows
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
    # print 'Max Point-to-Point Distance = %f' % np.sqrt(np.amax(sqDistMatrix))
    # print 'Mean Point-to-Point Distance = %f' % np.sqrt(np.mean(sqDistMatrix))
    # print np.shape(sqDistMatrix)
    #
    # distances = np.unique(sqDistMatrix)
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
    
    SDsquared = SD*SD
    # TwoSDsquared = 4*SDsquared

    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)
    
    cover_rows = [np.nonzero(t)[0] for t in C]

    return 0


def BuildModel(m):
    
    # DECLARE VARIABLES:
    # Facility Site binary decision variables X
    # Each has a coefficient of 1 in the objective
    sitesRange = range(numSites)
    X = m.addVars(sitesRange,
                  vtype=GRB.BINARY,
                  obj=np.ones(numSites),
                  name="X")
    
    # Define Coverage Constraints:
    for i in range(numDemands):
        m.addConstr(quicksum(X[j] for j in cover_rows[i]) >= 1)
    
    # The objective is to minimize the number of located facilities
    m.modelSense = GRB.MINIMIZE
    m.update()
    
    print 'Number of variables = %d' % m.numintvars
    print 'Number of constraints = %d' % m.numconstrs
    #m.printStats()
    
    print
    return 0

def SolveModel(m):
    """Solve the problem and print the solution."""
    m.Params.OutputFlag = 0
    m.Params.ResultFile = "output.sol"
    m.optimize()
    
    
def displaySolution(m, p, total_time):

    print 'Total problem solved in %f seconds' % total_time
    print
    # The objective value of the solution.
    print 'p = %d' % p
    print 'SD = %f' % SD
    # print the selected sites
    print
    j = 0    
    for v in m.getVars():
        if (v.x == 1.0):
            print "Site selected %s" % int(siteIDs[j])
        j += 1
    
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # plot solution 
    # plot.plotSolution(sites, X, range(numSites), SD)
    

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


def main(unused_argv):
    print ('---- LSCP with Gurobi -----')
    RunLSCP(SD)


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