"""Functions for kinetic Monte Carlo (Gillespie simulations)
"""
import numpy as np
import math
import random as rand

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
