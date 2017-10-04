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

def runProblem():
    #print sites
    #print sites[:,[0,3]]
    print sites[:,[1,2]]
    print np.shape(sites)
    p = 1
    computeCoverageMatrix(p, SD)
    
def computeCoverageMatrix(p, SD):
        
    #declare a couple variables
    global Nrows
    global Ncols
    global Nsize
    global facilityIDs
    global numSites
    
    # Convert Coordinates from Lat/Long to CONUS EqD Projection
    xyPointArray = sites[:,[1,2]]
    #A = [xyPointArray[i][:] for i in demandIDs]
    #B = [xyPointArray[j][:] for j in facilityIDs]
    A = xyPointArray
    B = A
    
    # Compute the distance matrix, using the squared distance
    sqDistMatrix = cdist(A, B,'sqeuclidean')
    print 'Max Point-to-Point Distance = %f' % np.sqrt(np.amax(sqDistMatrix))
    print 'Mean Point-to-Point Distance = %f' % np.sqrt(np.mean(sqDistMatrix))
    
#     if (p == -2):
#         p = 3
#     if (SD == -2):
#         SD = np.floor(np.sqrt(np.mean(sqDistMatrix))/13750.0)*1000
#
#     SDsquared = SD*SD
#     TwoSDsquared = 4*SDsquared
#
#     # Determine neighborhood of demands within SD of sites
#     C = (sqDistMatrix <= SDsquared).astype(int)
#
#     # Determine neighborhood of sites within SD of sites
#     if allFD3 == True:
#         SDist = C
#     else:
#         SDist = (cdist(B, B,'sqeuclidean') <= SDsquared).astype(int)
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
#     # Convert coverage to sparse matrix
#     Nrows,Ncols = np.nonzero(C.astype(bool))
#     Nsize = len(Nrows)
#
#     return [p, SD]
    return 0


# def computeCoverageMatrix(p, SD):
#
#     #declare a couple variables
#     global Nrows
#     global Ncols
#     global Nsize
#     global facilityIDs
#     global numSites
#
#     # Convert Coordinates from Lat/Long to CONUS EqD Projection
#     xyPointArray = GISOps.GetCONUSeqDprojCoords(js)
#     A = [xyPointArray[i][:] for i in demandIDs]
#     B = [xyPointArray[j][:] for j in facilityIDs]
#
#     # Compute the distance matrix, using the squared distance
#     sqDistMatrix = cdist(A, B,'sqeuclidean')
#     print 'Max Point-to-Point Distance = %f' % (np.sqrt(np.amax(sqDistMatrix))/1000.0)
#     print 'Mean Point-to-Point Distance = %f' % (np.sqrt(np.mean(sqDistMatrix))/1000.0)
#
#     if (p == -2):
#         p = 3
#     if (SD == -2):
#         SD = np.floor(np.sqrt(np.mean(sqDistMatrix))/13750.0)*1000
#
#     SDsquared = SD*SD
#     TwoSDsquared = 4*SDsquared
#
#     # Determine neighborhood of demands within SD of sites
#     C = (sqDistMatrix <= SDsquared).astype(int)
#
#     # Determine neighborhood of sites within SD of sites
#     if allFD3 == True:
#         SDist = C
#     else:
#         SDist = (cdist(B, B,'sqeuclidean') <= SDsquared).astype(int)
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
#     # Convert coverage to sparse matrix
#     Nrows,Ncols = np.nonzero(C.astype(bool))
#     Nsize = len(Nrows)
#
#     return [p, SD]
    

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


def main(unused_argv):
    runProblem()


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