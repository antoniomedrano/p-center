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

def Run_pCenter(p):
    
    """ Example of solving a 2-Center problem using brute force"""  
    
    start_time = time.time()
    
    distMatrix = computeDistanceMatrix()
    #print distMatrix
    
    SDmin, locations = SolveModel(p, distMatrix)
    
    total_time = time.time()-start_time
    #SDmin = m.objVal
    
    displaySolution(locations, p, SDmin, total_time)
    
    
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
    distMatrix = cdist(A, B,'euclidean')

    return distMatrix


def SolveModel(p, d):

    zbest = 1000000000000000.0

    for i in range(numSites):
        for j in range(i+1,numSites):
            for k in range(j+1,numSites):
                zmin = np.amax(np.array([d[i,:],d[j,:],d[k,:]]).min(0))
                # zmin = np.amax(np.amin([d[i,:],d[j,:],d[k,:]],axis=0)) # this is slower
                if (zmin < zbest):
                    zbest = zmin
                    s1, s2, s3 = i, j, k

    return zbest, [s1, s2, s3]
    
    
def displaySolution(locations, p, zbest, total_time):

    print 'Total problem solved in %f seconds' % total_time
    print
    # The objective value of the solution.
    print 'p = %d' % p
    print 'SD = %f' % zbest
    # print the selected sites
    print
    for j in locations:
        print "Site selected %s" % int(siteIDs[j])
    
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


def main(unused_argv):
    print ('---- P-Center with Gurobi -----')
    Run_pCenter(p)


""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 2 and len(sys.argv) <= 3:
    file = '../data/' + sys.argv[2]
    p = float(sys.argv[1])
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  elif len(sys.argv) > 1 and len(sys.argv) <= 2:
    p = float(sys.argv[1])
    file = '../data/swain.dat'
    print "Problem instance from: ", file
    read_problem(file)
    main(None)
  else:
    print "Please Pass: Service Distance; Data to Use"
    print "Problem not executed!"