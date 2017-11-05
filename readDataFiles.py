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

def readDat(file):
    """
    Read Rick Church's bespoke .dat file format
    It's essentially a CSV, with the following format for each line:
    ID, X, Y, population
    """
    
    i = 0
        
    # Use With Statement to automatically close the 'read file' when finished.
    with open(file,'r') as f:
        
        numSites =  sum(1 for _ in f)
        f.seek(0)
        sites = np.empty([numSites,4])        
        
        for line in f:       
            line = line.strip() # removes whitespace before and after the line
            
            # ignore empty lines
            if (len(line) == 0):
                sites = sites[:-1]
                numSites -= 1
                continue
            
            # ignore comments
            if (line[0] == '#' or line[0] == '%'):
                sites = sites[:-1]
                numSites -= 1
                continue            
                        
            row = line.split(",")
            # Set constraint coefficients
            sites[i,:] = row[0:4]
            i += 1
        
        return sites

