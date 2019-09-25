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
    
    plot.plotData(sites)
    
    print('%d locations' % numSites)
    print('Finished Reading File!')


def main(unused_argv):
    print('')
    #RunCBC_LSCPexampleCppStyleAPI(SD)
    #RunSCIP_LSCPexampleCppStyleAPI(SD)
    #RunBOP_LSCPexampleCppStyleAPI(SD)


""" Main will take in 3 arguments: p-Facilities; ServiceDistance; Data to Use  """
if __name__ == '__main__':
  if len(sys.argv) > 1:
    file = '../data/' + sys.argv[1]
    print("Problem instance from: ", file)
    read_problem(file)
    main(None)
  else:
    print("Problem not executed!")