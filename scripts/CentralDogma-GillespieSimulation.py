#!/usr/bin/env python

import os,sys
import numpy as np
import math
import random as rand
import matplotlib.pyplot as plt

def next_jump(total_propensity):
    xi = rand.uniform(1e-6,1.0)
    dt = -math.log(xi)/total_propensity

    return dt

def compute_propensities(central_dogma_rates,state):
    propensities = {}

    propensities['O_on'] = central_dogma_rates['k_on']*state['O_off']
    propensities['O_off'] = central_dogma_rates['k_off']*state['O_on']

    propensities['m'] = central_dogma_rates['k_m']*state['O_on']
    propensities['dm'] = central_dogma_rates['k_dm']*state['m']

    propensities['g'] = central_dogma_rates['k_g']*state['m']
    propensities['dg'] = central_dogma_rates['k_dg']*state['g']

    return propensities


def next_jump_and_event_type(propensities):
    total_wt = sum(list(propensities.values()))

    xi = rand.uniform(0,1.0)

    low_wt = 0.0

    for key in propensities.keys():
        if xi<=low_wt+(propensities[key]/total_wt) and xi>low_wt:
            event = key
            break
        else:
            low_wt += propensities[key]/total_wt

    xi = rand.uniform(1e-8,1.0)

    dt = -math.log(xi)/total_wt

    event_prob = propensities[event]/total_wt

    return event, dt, event_prob

def update_state(event,state):
    if event=='O_on':
        state['O_on'] = 1
        state['O_off'] = 0
    elif event=='O_off':
        state['O_on'] = 0
        state['O_off'] = 1
    elif event=='m':
        state['m'] += 1
    elif event=='dm':
        state['m'] += -1
    elif event=='g':
        state['g'] += 1
    else:
        state['g'] += -1

    return state

def reverse_event_prob(event,state,central_dogma_rates):
    rev_propensities = compute_propensities(central_dogma_rates,state)
    total_wt = sum(list(rev_propensities.values()))

    if event=='O_on':
        prob = rev_propensities['O_off']/total_wt
    elif event=='O_off':
        prob = rev_propensities['O_on']/total_wt
    elif event=='m':
        prob = rev_propensities['dm']/total_wt
    elif event=='dm':
        prob = rev_propensities['m']/total_wt
    elif event=='g':
        prob = rev_propensities['dg']/total_wt
    else:
        prob = rev_propensities['g']/total_wt

    return prob

data_directory = '/Users/sns9/Research/Translation/tpower16/'

try:
    os.chdir(data_directory)
except OSError:
    os.mkdir(data_directory)
    os.chdir(data_directory)

# leaky expression coefficient
l = 0.01
s = 1.0

# Input code length
lenX = 4
X = np.linspace(0,1.0,int(2**lenX))
#pr_array = np.zeros(shape=(X.shape[0],3))

# Number of samples
n_samples = 100000

# Frequency factor for k_on and k_off
alpha = 1.0

# Transcription and translation parameters
k_m = 0.5
k_dm = 0.5

translation_power = 1.0
translation_power_tag = str(translation_power).replace('.','-')
# Protein decay rates
k_dgs = [0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002]

expression_directory = translation_power_tag+'_samples'

try:
    os.chdir(expression_directory)
except OSError:
    os.mkdir(expression_directory)
    os.chdir(expression_directory)

for k_dg in k_dgs:
    k_g = translation_power*k_dg

    outstring = str(k_m)+','+str(k_dm)+','+str(k_g)+','+str(k_dg)+','+str(k_dm/k_dg)+','

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
    dissipation = {}

    g_sample_array = np.zeros(shape=(n_samples,X.shape[0]))

    for i in range(0,X.shape[0]):
        current_time = 0.0
        event_counter = 0
        last_sample_event = 0

        # Transcription on or off
        k_on = alpha*((1-l)*X[i] + l)
        k_off = alpha*(1 - X[i])*(1 - l)

        # Initial transcript count
        m = int((k_m/k_dm)*k_on/(k_on + k_off))

        #m_samples[i] = [m]

        g = int((k_g/k_dg)*(k_m/k_dm)*k_on/(k_on + k_off))

        #g_samples[i] = [g]
        g_sample_array[0,i] = state['g']

        current_time = 0
        samples = 1

        state['O_on'] = 0
        state['O_off'] = 1
        state['m'] = m
        state['g'] = g

        central_dogma_rates['k_on'] = k_on
        central_dogma_rates['k_off'] = k_off

        propensities = compute_propensities(central_dogma_rates,state)

        t_sample = 10.0
        next_sample_time = t_sample

        dissipation[i] = 0.0

        while samples<n_samples:
            propensities = compute_propensities(central_dogma_rates,state)

            event, dt, event_prob = next_jump_and_event_type(propensities)

            state = update_state(event,state)

            rev_prob = reverse_event_prob(event,state,central_dogma_rates)

            if event_prob>0.0 and rev_prob>0.0:
                dissipation[i] += math.log(event_prob/rev_prob)

            current_time += dt
            event_counter += 1

            if current_time>=next_sample_time:
                next_sample_time = current_time + t_sample - (current_time - next_sample_time)

                g_sample_array[samples,i] = state['g']

                samples += 1
                last_sample_event = event_counter

        dissipation[i] *= 1.0/current_time

        print(X[i],' is completed: ',current_time,event_counter)

    tfile = 'T'+str(k_dm/k_dg).replace('.','_')+'.csv'

    np.savetxt(tfile,g_sample_array,delimiter=',')
