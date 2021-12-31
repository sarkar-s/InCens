Computation of protein-level information gain curves
====================================================

Computes the channel capacity for the protein expression level from the
Gillespie simulation results and further corrects for overesimation due
to finite number of samples.

.. code:: ipython3

    import os,sys
    import numpy as np
    import scipy.stats as st
    import scripts.BA_C as bao
    import math
    import random as rand
    import matplotlib.pyplot as plt
    import glob
    import pandas as pd
    import scipy.stats as st
    from matplotlib import rcParams
    import matplotlib.gridspec as gridspec

.. code:: ipython3

    def get_subsample_data(this_data,df):
        sub_data = np.zeros(shape=(int(this_data.shape[0]/df),this_data.shape[1]))
        
        for i in range(0,this_data.shape[1]):
            sub_data[:,i] = rand.choices(this_data[:,i],k=sub_data.shape[0])
            
        return sub_data

Select fractions for which to sample (with replacements) the stochastic
simulation data and the number of replicates for each fraction.

.. code:: ipython3

    data_fractions = [1,2,5]
    reps = 10

Read the protein expression levels from the Gillespie simulation of the
central dogma reaction system.

.. code:: ipython3

    data_directory = './simulation_results/test_samples'
    os.chdir(data_directory)
    
    c_files = glob.glob('T*.csv')

.. code:: ipython3

    t_values = []
    datas = {}
    diss_datas = {}
    
    for f in c_files:
        t = f.replace('T','').replace('.csv','').replace('_','.')
    
        datas[t] = pd.read_csv(f,header=None).to_numpy()
        
        t_values.append(float(t))
    
    t_values.sort()
    
    t_and_c = np.zeros(shape=(len(t_values),2))
    t_and_c[:,0] = np.array(t_values)

Calculate the channel capacity for all data fractions and replicates and
then linearly fit the channel capacity values vs the data fractions. The
y-intercept of the linear fit is the channel capacity after correcting
for the finite-sampling bias.

.. code:: ipython3

    max_bins = 32
    
    cs = []
    
    print('T \t c(X;g)')
    
    for i in range(0,t_and_c.shape[0]):
        d = str(t_and_c[i,0])
        
        gmin, gmax = np.min(datas[d]),np.max(datas[d])
        
        if t_and_c[i,0]<10.0:
            bins = max_bins/8
        elif t_and_c[i,0]<100.0:
            bins = max_bins/4
        elif t_and_c[i,0]<1000.0:
            bins = max_bins/2
     
        bin_size = int(min(bins,int(gmax - gmin + 1)))
                
        g_edges = np.linspace(gmin,gmax,bin_size+1)
        
        this_c_set = np.zeros(shape=(len(data_fractions)*reps,2))
        
        k = 0
        
        for df in data_fractions:
            for r in range(0,reps):
                sub_data = get_subsample_data(datas[d],df)
                
                g_pdfs = np.zeros(shape=(datas[d].shape[1],bin_size))
    
                for j in range(0,datas[d].shape[1]):
                    ghist, bin_edges = np.histogram(sub_data[:,j],bins=g_edges)
    
                    g_pdfs[j,:] = ghist/np.sum(ghist)
        
                c, e, p = bao.get_CC(g_pdfs)
                
                this_c_set[k,0] = df
                this_c_set[k,1] = c
                
                k += 1
                
        res = st.linregress(this_c_set[:,0],this_c_set[:,1])
        
        c = float("{:.2f}".format(res.intercept))
        
        t_and_c[i,1] = c
        
        print(d,'\t',c)
        
    sorted_t_and_c = t_and_c[tuple([np.argsort(t_and_c[:,0])])]

Plot and save the information gain curve.

.. code:: ipython3

    outfile = 'C-summary.csv'
    
    np.savetxt(outfile,sorted_t_and_c,delimiter=',',header='T,c_g',comments='')
    
    fig = plt.figure(figsize=(5,4))
    
    plt.plot(t_and_c[:,0],t_and_c[:,1],marker='.',ms=10)
    plt.xscale('log')
    plt.xlabel(r'Integration time, $T=k_{d,m}/k_{d,g}$',size=16)
    plt.ylabel(r'$c(X;g)$ (bits)',size=16)
    plt.tick_params(labelsize=16)

