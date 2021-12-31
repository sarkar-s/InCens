"""
Stochastic simulator for central dogma reaction system.
"""

import os,sys
import numpy as np
import scipy.stats as st
import math
import random as rand
from optparse import OptionParser

#import kMC as kMC

from simulation_data import simulation_data


"""
Receives the central dogma rate constants data name and number of samples to capture from the stochastic protein expression trajectory.
Both parameters are recieved as command line arguments.

Args:
    dataname (string): Data identifier within quotes, e.g. "Homo sapiens", which is then read from the simulation_data.py module.

    samples (int): Number of protein expression values to sample from the stochastic expression trajectory.
"""

parser = OptionParser()
parser.add_option("-d","--data",dest="dataname",help="data field (in quotes) containing the simulation parameters")
parser.add_option("-n","--samples",dest="samples",help="number of protein expression values to sample during the simulation")

(options, args) = parser.parse_args()

selected_data = options.dataname

# Number of samples
n_samples = int(options.samples)

data_directory = '../simulation_results/'

try:
    os.chdir(data_directory)
except OSError:
    os.mkdir(data_directory)
    os.chdir(data_directory)

expression_directory = selected_data.replace(' ','-')+'_samples'

try:
    os.chdir(expression_directory)
except OSError:
    os.mkdir(expression_directory)
    os.chdir(expression_directory)


# Frequency factor for k_on and k_off
alpha = 1.0
# leakiness
l = 0.01

# Transcription parameters
k_m = simulation_data[selected_data]['k_m']
k_dm = simulation_data[selected_data]['k_dm']
k_dgs = simulation_data[selected_data]['k_dg']
translation_power = simulation_data[selected_data]['translation power']

for k_dg_i in range(0,len(k_dgs)):
    k_dg = k_dgs[k_dg_i]
    k_g = translation_power*k_dg

    state = {}
    state['O_on'] = 0
    state['O_off'] = 1
    state['m'] = 0
    state['g'] = 0

    central_dogma_rates = {}
    central_dogma_rates['k_m'] = k_m
    central_dogma_rates['k_dm'] = k_dm
    central_dogma_rates['k_g'] = k_g
    central_dogma_rates['k_dg'] = k_dg

    g_samples = {}
    m_samples = {}

    print(simulation_data[selected_data]['input size'][k_dg_i])

    X = np.linspace(0,1.0,int(simulation_data[selected_data]['input size'][k_dg_i]))

    g_sample_array = np.zeros(shape=(n_samples,X.shape[0]))
    dissipation = np.zeros(shape=(X.shape[0],))

    t_sample = int(simulation_data[selected_data]['sample interval'][k_dg_i])

    for i in range(0,X.shape[0]):
        current_time = 0.0
        event_counter = 0
        last_sample_event = 0

        # Transcription on or off
        k_on = alpha*((1-l)*X[i] + l)
        k_off = alpha*(1 - X[i])*(1 - l)

        # Initial transcript count
        m = int((k_m/k_dm)*k_on/(k_on + k_off))

        m_samples[i] = [m]

        g = int((k_g/k_dg)*(k_m/k_dm)*k_on/(k_on + k_off))

        g_samples[i] = [g]
        g_sample_array[0,i] = state['g']

        state['O_on'] = 0
        state['O_off'] = 1
        state['m'] = m
        state['g'] = g
        g_sample_array[0,i] = state['g']

        current_time = 0
        samples = 1

        central_dogma_rates['k_on'] = k_on
        central_dogma_rates['k_off'] = k_off

        propensities = kMC.compute_propensities(central_dogma_rates,state)

        next_sample_time = t_sample

        dissipation[i] = 0.0

        while samples<n_samples:
            propensities = kMC.compute_propensities(central_dogma_rates,state)

            event, dt, event_prob = kMC.next_jump_and_event_type(propensities)

            state = kMC.update_state(event,state)

            rev_prob = kMC.reverse_event_prob(event,state,central_dogma_rates)

            if event_prob>0.0 and rev_prob>0.0:
                dissipation[i] += math.log(event_prob/rev_prob)

            current_time += dt
            event_counter += 1

            if current_time>=next_sample_time:
                m_samples[i].append(state['m'])
                g_samples[i].append(state['g'])

                next_sample_time = current_time + t_sample - (current_time - next_sample_time)

                g_sample_array[samples,i] = state['g']

                samples += 1
                last_sample_event = event_counter

        dissipation[i] *= 1.0/current_time

        print(X[i],' is completed: ',current_time,event_counter)

    Tvalue = float("{:.1f}".format(k_dm/k_dg))

    tfile = 'T'+str(Tvalue).replace('.','_')+'.csv'

    np.savetxt(tfile,g_sample_array,delimiter=',')
