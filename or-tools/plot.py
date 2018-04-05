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
    
    mSize = 4
    
    plt.figure(figsize=(6,6))
    plt.plot(sites[:,1], sites[:,2], 'ko', markersize = mSize)
    plt.axis('equal')
    plt.show()

def plotSolution(sites, X, cols, SD):
    
    mSize = 4
    
    plt.figure(figsize=(6,6))
    plt.plot(sites[:,1], sites[:,2], 'ko', markersize = mSize)
    numSites = len(X)
    for j in range(numSites):
        if (X[j].SolutionValue() == 1.0):
            plt.plot(sites[cols[j],1], sites[cols[j],2], 'bo', markersize = mSize)
            circle = plt.Circle((sites[cols[j],1], sites[cols[j],2]), SD, color='b', fill=False)
            plt.gcf().gca().add_artist(circle)
    
    plt.axis('equal')
    plt.show()
    
def plotSolutionE(sites, essential, X, cols, SD):
    
    mSize = 4
    
    plt.figure(figsize=(6,6))
    plt.plot(sites[:,1], sites[:,2], 'ko', markersize = mSize)
    numSites = len(X)
    for j in range(numSites):
        if (X[j].SolutionValue() == 1.0):
            plt.plot(sites[cols[j],1], sites[cols[j],2], 'bo', markersize = mSize)
            circle = plt.Circle((sites[cols[j],1], sites[cols[j],2]), SD, color='b', fill=False)
            plt.gcf().gca().add_artist(circle)
    for j in essential:
        plt.plot(sites[j,1], sites[j,2], 'ro', markersize = mSize)
        circle = plt.Circle((sites[j,1], sites[j,2]), SD, color='b', fill=False)
        plt.gcf().gca().add_artist(circle)
    
    plt.axis('equal')
    plt.show()
    
def plotTradeoff(file, solution):

    mSize = 2
    rows, cols = solution.shape

    plt.figure(figsize=(10,6))
    plt.step(solution[:,1], solution[:,0], where='pre')
    plt.plot(solution[:,1], solution[:,0], 'ko', markersize = mSize)
    plt.ylim(0, solution[rows-1,0]+1)
    plt.xlim(0, solution[0,1]+1)
    #plt.xscale('log')
    plt.ylabel('p-facilities')
    plt.xlabel('Service Distance')    
    plt.title('Complete P-Center Trade-Off Frontier, ' + file)
    plt.show()
    
    # plt.figure(figsize=(10,6))
    # plt.step(solution[:,0], solution[:,1], where='post')
    # plt.plot(solution[:,0], solution[:,1], 'ko', markersize = mSize)
    # plt.xlim(0, solution[rows-1,0]+1)
    # plt.ylim(0, solution[0,1]+1)
    # #plt.yscale('log')
    # plt.xlabel('p-facilities')
    # plt.ylabel('Service Distance')
    # plt.title('Complete P-Center Trade-Off Frontier')
    # #plt.axis('equal')
    # plt.show()

