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
    plt.plot(sites[:,1], sites[:,2], 'bo')
    plt.axis('equal')
    plt.show()
    
def plotSolution(sites, X, SD):
    
    plt.figure(figsize=(6,6))
    plt.plot(sites[:,1], sites[:,2], 'bo')
    plt.axis('equal')
    plt.show()