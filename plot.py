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
import matplotlib.pyplot as plt
import numpy as np

def plotData(sites):
    
    plt.figure(figsize=(6,6))
    plt.plot(sites[:,1], sites[:,2], 'ko')
    plt.axis('equal')
    plt.show()

def plotSolution(sites, X, cols, SD):
    
    plt.figure(figsize=(6,6))
    plt.plot(sites[:,1], sites[:,2], 'ko')
    numSites = len(X)
    for j in range(numSites):
        if (X[j].SolutionValue() == 1.0):
            plt.plot(sites[cols[j],1], sites[cols[j],2], 'bo')
            circle = plt.Circle((sites[cols[j],1], sites[cols[j],2]), SD, color='b', fill=False)
            plt.gcf().gca().add_artist(circle)
    
    plt.axis('equal')
    plt.show()
    
def plotSolutionE(sites, essential, X, cols, SD):
    
    plt.figure(figsize=(6,6))
    plt.plot(sites[:,1], sites[:,2], 'ko')
    numSites = len(X)
    for j in range(numSites):
        if (X[j].SolutionValue() == 1.0):
            plt.plot(sites[cols[j],1], sites[cols[j],2], 'bo')
            circle = plt.Circle((sites[cols[j],1], sites[cols[j],2]), SD, color='b', fill=False)
            plt.gcf().gca().add_artist(circle)
    for j in essential:
        plt.plot(sites[j,1], sites[j,2], 'ro')
        circle = plt.Circle((sites[j,1], sites[j,2]), SD, color='b', fill=False)
        plt.gcf().gca().add_artist(circle)
    
    plt.axis('equal')
    plt.show()