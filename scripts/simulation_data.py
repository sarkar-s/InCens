#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os,sys
import numpy as np
import scipy.stats as st
import math
import random as rand
import matplotlib.pyplot as plt

simulation_data = {}

"""All rate constants are in min^{-1}
"""

# Data for the 3 model organisms and homo sapiens

simulation_data['Escherichia coli'] = {}
simulation_data['Saccharomyces cerevisiae'] = {}
simulation_data['Mus musculus'] = {}
simulation_data['Homo sapiens'] = {}

simulation_data['Escherichia coli']['k_m'] = 0.02
simulation_data['Escherichia coli']['k_dm'] = 0.15
simulation_data['Escherichia coli']['k_dg'] = [0.15,0.075,0.03,0.015,0.0075,0.003,0.0015,0.00075,0.0003,0.00015]
simulation_data['Escherichia coli']['translation power'] = 600
simulation_data['Escherichia coli']['input size'] = [16,16,16,16,32,32,32,32,64,64]
simulation_data['Escherichia coli']['sample interval'] = [10,10,10,10,10,10,10,10,10,10]

simulation_data['Saccharomyces cerevisiae']['k_m'] = 0.31
simulation_data['Saccharomyces cerevisiae']['k_dm'] = 0.1
simulation_data['Saccharomyces cerevisiae']['k_dg'] = [0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001]
simulation_data['Saccharomyces cerevisiae']['translation power'] = 400
simulation_data['Saccharomyces cerevisiae']['input size'] = [16,16,16,16,32,32,32,32,64,64]
simulation_data['Saccharomyces cerevisiae']['sample interval'] = [10,10,10,10,10,10,10,10,10,10]

simulation_data['Mus musculus']['k_m'] = 0.02
simulation_data['Mus musculus']['k_dm'] = 0.002
simulation_data['Mus musculus']['k_dg'] = [0.002,0.001,0.0004,0.0002,0.0001,0.00004,0.00002,0.00001,0.000004,0.000002]
simulation_data['Mus musculus']['translation power'] = 25000
simulation_data['Mus musculus']['input size'] = [16,16,16,16,32,32,64,64,128,128]
simulation_data['Mus musculus']['sample interval'] = [10,10,10,10,10,10,100,100,200,200]

simulation_data['Homo sapiens']['k_m'] = 0.01
simulation_data['Homo sapiens']['k_dm'] = 0.002
simulation_data['Homo sapiens']['k_dg'] = [0.002,0.001,0.0005,0.0002,0.0001,0.00005,0.00002,0.00001,0.000005,0.000002]
simulation_data['Homo sapiens']['translation power'] = 20000
simulation_data['Homo sapiens']['input size'] = [16,16,16,16,32,32,64,64,128,128]
simulation_data['Homo sapiens']['sample interval'] = [10,10,10,10,10,10,100,100,200,200]

# Data for the protein-level information gain curves in Fig. 2(b)
simulation_data['test'] = {}

simulation_data['test']['k_m'] = 0.5
simulation_data['test']['k_dm'] = 0.5
simulation_data['test']['k_dg'] = [0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002]
simulation_data['test']['translation power'] = 1 # Translation power set is 1,10,100,1000,10000
simulation_data['test']['input size'] = [16,16,16,16,16,64,64,64,64,64,64]
simulation_data['test']['sample interval'] = [10,10,10,10,10,10,10,10,10,10,10]
