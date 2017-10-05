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
import numpy as np
from scipy.sparse import csc_matrix
from scipy.spatial.distance import cdist
from ortools.linear_solver import pywraplp

def RunLSCPCppStyleAPI(optimization_problem_type, p, SD):
    
    """ Example of simple MCLP program with the C++ style API."""
    solver = pywraplp.Solver('RunIntegerExampleCppStyleAPI', optimization_problem_type)
    
    
    #print sites
    #print np.shape(sites)
    computeCoverageMatrix(p, SD)
    
    
def computeCoverageMatrix(p, SD):
        
    #declare a couple variables
    global distances
    global numDemands
    global numSites
    global Nrows
    global Ncols
    global Nsize
    
    # for now, all demands are also sites
    allFD3 = True
    facilityIDs = range(len(sites))
    
    # Convert Coordinates from Lat/Long to CONUS EqD Projection
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
    
    print colmax
    print minmax
    
    print np.where(distances==minmax)
    
    SDsquared = SD*SD
    TwoSDsquared = 4*SDsquared

    # Determine neighborhood of demands within SD of sites
    C = (sqDistMatrix <= SDsquared).astype(int)

    # Determine neighborhood of sites within SD of sites
    if allFD3 == True:
        SDist = C
    else:
        SDist = (cdist(B, B,'sqeuclidean') <= SDsquared).astype(int)
#
#     start_time = time.time()
#     C, columns = dominationTrim(C, SDist)
#     print 'Domination time = %f' % (time.time()-start_time)
#
#     # shorten the facility data sets
#     cols = np.nonzero(columns)[0]
#     facilityIDs = [facilityIDs[j] for j in cols]
#     numSites = len(facilityIDs)
#
    # Convert coverage to sparse matrix
    Nrows,Ncols = np.nonzero(C.astype(bool))
    Nsize = len(Nrows)
#
#     return [p, SD]
    return 0

def BuildModel(solver, X, Y, p):
    
    infinity = solver.infinity()
    
    # DECLARE CONSTRAINTS:
    # declare demand coverage constraints (binary integer: 1 if UNCOVERED, 0 if COVERED)
    c1 = [None]*numDemands
    
    
    # declare the objective
    objective = solver.Objective()
    objective.SetMinimization()
    
    # Add potential facility sites to the p constraint
    for j in range(numDemands):
        # initialize the X variables as Binary Integer (Boolean) variables
        name = "X,%d" % facilityIDs[j]
        X[j] = solver.BoolVar(name)
        c2.SetCoefficient(X[j],1)
    
    # if facility is fixed into the solution, add a constraint to make it so
    for k in range(numForced):
          c3[k] = solver.Constraint(1,1)
          c3[k].SetCoefficient(X[forcedFacilities[k]],1)
    
    # add demands to the objective and coverage constraints
    for i in range(numDemands):
        name = "Y,%d" % demandIDs[i]
        # initialize the Y variables as Binary Integer (Boolean) variables
        Y[i] = solver.BoolVar(name)
        # Set the Objective Coefficients for the population * Demand Variable (Yi)
        objective.SetCoefficient(Y[i],demandPop[i])
        # Covering constraints
        c1[i] = solver.Constraint(1, solver.infinity())
        c1[i].SetCoefficient(Y[i],1)

    # add facility coverages to the coverage constraints
    for k in range(Nsize):
        c1[Nrows[k]].SetCoefficient(X[Ncols[k]],1)
    
    print 'Number of variables = %d' % solver.NumVariables()
    print 'Number of constraints = %d' % solver.NumConstraints()
    print
    return 0
    

def read_problem(file):
    global numSites
    global sites
    
    print 'readFile({0})'.format(file)
    
    lineCount = 0
    i = 0
    
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

def RunSCIPPMedianExampleCppStyleAPI(p, SD):
    if hasattr(pywraplp.Solver, 'SCIP_MIXED_INTEGER_PROGRAMMING'):
        Announce('SCIP', 'C++ style API')
        RunLSCPCppStyleAPI(pywraplp.Solver.SCIP_MIXED_INTEGER_PROGRAMMING, p, SD)

def RunCBCMCLPexampleCppStyleAPI(p, SD):
    if hasattr(pywraplp.Solver, 'CBC_MIXED_INTEGER_PROGRAMMING'):
        Announce('CBC', 'C++ style API')
        RunLSCPCppStyleAPI(pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING, p, SD)

def RunBOPMCLPexampleCppStyleAPI(p, SD):
    if hasattr(pywraplp.Solver, 'BOP_INTEGER_PROGRAMMING'):
        Announce('BOP', 'C++ style API')
        RunLSCPCppStyleAPI(pywraplp.Solver.BOP_INTEGER_PROGRAMMING, p, SD)


def main(unused_argv):
    p = 1
    RunCBCMCLPexampleCppStyleAPI(p, SD)


""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 2 and len(sys.argv) <= 3:
    file = sys.argv[2]
    SD = float(sys.argv[1])
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  elif len(sys.argv) > 1 and len(sys.argv) <= 2:
    SD = long(sys.argv[1])
    file = r'./data/swain.txt'
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  else:
    print "Please Pass: Service Distance; Data to Use"
    print "Problem not executed!"