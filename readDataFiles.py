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
import re
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


def readTSP(file):
    """
    Read a file from the TSP library .tsp file format
    """
    
    data = read_tsp_data(file)
    numSites = int(detect_dimension(data))
    sites = np.ones([numSites,4])        
    i = 0
    dataPortion = False
    for line in data:
        if dataPortion == True:
            row = line.split()
            # Set constraint coefficients
            #print row
            sites[i,0:3] = row[0:3]
            i += 1
            if i == numSites:
                break
        else:
            if line.startswith('NODE_COORD_SECTION'):
                dataPortion = True
                continue
    return sites
    
    
def read_tsp_data(tsp_name):
	tsp_name = tsp_name
	with open(tsp_name) as f:
		content = f.read().splitlines()
		cleaned = [x.lstrip() for x in content if x != ""]
		return cleaned
"""
We return a list like 
['NAME: ulysses16.tsp',
'TYPE: TSP',
'COMMENT: Odyssey of Ulysses (Groetschel/Padberg)',
'DIMENSION: 16',
'EDGE_WEIGHT_TYPE: GEO',
'DISPLAY_DATA_TYPE: COORD_DISPLAY',
'NODE_COORD_SECTION',
'1 38.24 20.42',
'2 39.57 26.15',
'3 40.56 25.32',
................
'EOF']
"""

"""
Check for the line DIMENSION in the file and keeps the numeric value
"""
def detect_dimension(in_list):
	non_numeric = re.compile(r'[^\d]+')
	for element in in_list:
		if element.startswith("DIMENSION"):
			return non_numeric.sub("",element)
