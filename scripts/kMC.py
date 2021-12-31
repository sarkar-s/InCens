"""
Functions for kinetic Monte Carlo simulation (Gillespie algorithm)
"""
import numpy as np
import math
import random as rand

def compute_propensities(central_dogma_rates,state):
    r"""Determines the propensities of the central dogma reaction events as a function of the central dogma system state.
    The states is defined using :math:`O`- operator state, :math:'m'- transcript count, and :math:`g` protein copy number. The propensities
    of the following events are computed: operator switch from Off to On, operator switch from On to Off, transcription,
    transcript decay, translation, and protein decay.

    Parameters
    ----------
    central_dogma_rates: dict
        Dictionary containing the rate constants, :math:`k_{on},k_{off},k_m,k_{dm},k_g,` and :math:`k_{dg}`.

    state: dict
        Dictionary containing the operator state value, the transcript count, and the protein count.

    Returns
    -------
    propensities: dict
        Propensities of all the central dogma reactions.
    """

    propensities = {}

    propensities['O_on'] = central_dogma_rates['k_on']*state['O_off']
    propensities['O_off'] = central_dogma_rates['k_off']*state['O_on']

    propensities['m'] = central_dogma_rates['k_m']*state['O_on']
    propensities['dm'] = central_dogma_rates['k_dm']*state['m']

    propensities['g'] = central_dogma_rates['k_g']*state['m']
    propensities['dg'] = central_dogma_rates['k_dg']*state['g']

    return propensities

def next_jump_and_event_type(propensities):
    r"""
    Determines the time interval till next reaction using the total propensity of all the reaction events and the reaction event.

    Parameters
    ----------
    propensities: dict
        Propensities all the reactions in the model system.

    Returns
    -------
    event_type: string
        Name of the reaction event that will occur in the next step of stochastic simulation.

    dt: float
        Time interval to the next reaction event.

    event_prob: float
        Probability of occurence of the selected reaction event.
    """

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
    r"""Updates the state of the central dogma reaction system based on the selected reaction event.

    Parameters
    ----------
    event: string
        Propensities all the reactions in the model system.

    state: dict
        Dictionary containing the operator state value, and the transcript and protein counts.

    Returns
    -------
    state: dict
        Dictionary containing the updated operator state value, and the transcript and protein counts,
        based on the selected event.
    """

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
    r"""Determines the propbability of the occurence of the reverse of the selected event.
    If the event probability of the selected event is :math:`P(S_2|S_1)`, then this function calculates :math:`P(S_1|S_2)`.

    Parameters
    ----------
    event: string
        Propensities all the reactions in the model system.

    state: dict
        Dictionary containing the operator state value, and the transcript and protein counts.

    central_dogma_rates: dict
        Dictionary containing the rate constants, :math:`k_{on},k_{off},k_m,k_{dm},k_g,` and :math:`k_{dg}`.

    Returns
    -------
    prob: float
        Probability of occurence of the reverse of the selected reaction event.
    """

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
