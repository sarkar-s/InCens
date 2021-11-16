#!/usr/bin/env python
# coding: utf-8

import os,sys
import numpy as np
import scipy.stats as st
import math
import random as rand
from BA_C import BA
from optparse import OptionParser

from simulation_data import simulation_data

parser = OptionParser()
parser.add_option("-d","--data",dest="dataname",help="data field (in quotes) containing the simulation parameters")

(options, args) = parser.parse_args()

selected_data = options.dataname

data_directory = '../simulation_results/'

try:
    os.chdir(data_directory)
except OSError:
    os.mkdir(data_directory)
    os.chdir(data_directory)

c_file = selected_data.replace(' ','-')+'-C_ideal.csv'

# Frequency factor for k_on and k_off
alpha = 1.0
# leakiness
l = 0.01

# Transcription parameters
k_m = simulation_data[selected_data]['k_m']
k_dm = simulation_data[selected_data]['k_dm']
k_dgs = simulation_data[selected_data]['k_dg']
translation_power = simulation_data[selected_data]['translation power']
all_cs = np.zeros(shape=(len(k_dgs),2))

# Blahut-Arimoto algorithm to compute channel capacity
bao = BA()

for k_dg_i in range(0,len(k_dgs)):
    k_dg = k_dgs[k_dg_i]
    k_g = translation_power*k_dg

    all_cs[k_dg_i,0] = k_dm/k_dg
    t = all_cs[k_dg_i,0]

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

    X = np.linspace(0,1.0,int(simulation_data[selected_data]['input size'][k_dg_i]))

    # Transcript expression distribution
    r_params = np.zeros(shape=(X.shape[0],2))

    for i in range(0,X.shape[0]):
        k_on = alpha*((1-l)*X[i] + l)
        k_off = alpha*(1 - X[i])*(1 - l)

        m = (k_m/k_dm)*k_on/(k_on + k_off)

        b = 1 + (k_dm*k_off*m)/(k_on*(k_on+k_off+k_dm))

        beta = b - 1

        if beta>0.0:
            r_params[i,0] = m/beta
            r_params[i,1] = (b - 1.0)/b
        else:
            r_params[i,0] = m
            r_params[i,1] = 0.0

    # create bins
    if r_params[0,1]>0.0:
        g_min = st.nbinom.ppf(0.001,r_params[0,0]*t,1-r_params[0,1])
    else:
        g_min = st.poisson.ppf(0.001,r_params[0,0]*t)

    if r_params[-1,1]>0.0:
        g_max = st.nbinom.ppf(0.999,r_params[-1,0]*t,1-r_params[-1,1])
    else:
        g_max = st.poisson.ppf(0.999,r_params[-1,0]*t)

    bin_size = min(256,g_max - g_min + 1)

    g_bin_edges = np.linspace(g_min,g_max+1,bin_size+1)
    g_locs_i = g_bin_edges.astype(int)
    g_pdfs = np.zeros(shape=(r_params.shape[0],g_locs_i.shape[0]))

    for i in range(0,r_params.shape[0]):
        if r_params[i,1]>0.0:
            p = r_params[i,1]
            r = r_params[i,0]*t

            l_cdf = st.nbinom.cdf(g_locs_i,r,1-p,0)
            l_pdf = l_cdf
            l_pdf[1:] = l_pdf[1:] - l_pdf[:-1]
        else:
            r = r_params[i,0]*t
            l_cdf = st.poisson.cdf(g_locs_i,r,0)
            l_pdf = l_cdf
            l_pdf[1:] = l_pdf[1:] - l_pdf[:-1]

        g_pdfs[i,:] = l_pdf/np.sum(l_pdf)

    bao.set_response(g_pdfs)
    c_g, e, p = bao.get_CC()

    all_cs[k_dg_i,1] = c_g

    np.savetxt(c_file,all_cs,delimiter=',')
