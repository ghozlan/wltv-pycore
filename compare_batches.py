# -*- coding: utf-8 -*-
"""
Created on Fri May 08 20:30:11 2015

@author: Hassan
"""

#%%
# get MATLAB results
sim_filename_pattern = '[0-9]{6}_master\.mat$'
sim_filename_pattern = '[0-9]{6}_marx\.mat$'
sim_filename_pattern = '[0-9]{6}_py$'
sim_filename_pattern = 'thinice_py$'
sim_filename_pattern = 'sarajevo_py$'
sim_filename_pattern = '\.mat$'
import load_batch
RESULTS1 = load_batch.load_batch(sim_filename_pattern = '\.mat$', file_type='mat')
RESULTS2 = load_batch.load_batch(sim_filename_pattern = '_py$')
#RESULTS2 = load_batch.load_batch(sim_filename_pattern = 'marx.mat$', file_type='mat')

#%%
from matplotlib.pyplot import *
import numpy

counter = 0
ERR = {}

#CH_LABELS = ['A', 'B', 'C', 'D', 'E']
#SCH_LABELS = [1, 2, 3]
CH_LABELS = ['A', 'E']
SCH_LABELS = [1, 2]
from itertools import product # Caterzian/Cartesian product
for ch_index, sch_index in product(CH_LABELS, SCH_LABELS):
    #print('Channel %s, Scheme %d' %((str(ch_index), sch_index)))
    
    for RX in ['OPT','EB','SB_JLD','SB_ILD']:
        r1 = RESULTS1[('CHANNEL',ch_index)][('SCHEME',sch_index)][RX]
        r2  = RESULTS2 [('CHANNEL',ch_index)][('SCHEME',sch_index)][RX]
        r_error = numpy.abs(r1 - r2)
        #ERR[counter] = mean(r_error)
        ERR[(ch_index, sch_index, RX)] = numpy.mean(r_error)
        counter += 1
        #plot(SNRdB_vec,r_error)

#for i in ERR: 
#    if ERR[i]>1e-1: print ERR[i], '->', i

#%%
RX_LIST = ['OPT','EB','SB_JLD','SB_ILD']
sim_list = [(item[0][0],item[0][1],item[1]) for item in product(product(CH_LABELS, SCH_LABELS), RX_LIST)]
#print sim_list
  
for sim_idx in sim_list:
    if ERR[sim_idx]>1e-2:
        print '%.5f' % ERR[sim_idx], '->', sim_idx      

import sys
#sys.exit()
#%%

for sim_idx in sim_list: print ERR[sim_idx], '->', sim_idx
#sys.exit()
#%%
SNRdB_vec = RESULTS1['SNRdB']
i = 0
for sim in sim_list:
    ch_index,sch_index,RX = sim
    #ch_index,sch_index,RX = the_list[sim_idx]
    #ch_index,sch_index,RX = ('A', 2, 'SB_ILD')
    r1 = RESULTS1[('CHANNEL',ch_index)][('SCHEME',sch_index)][RX]
    r2  = RESULTS2 [('CHANNEL',ch_index)][('SCHEME',sch_index)][RX]
    i += 1
    figure(i)
    plot(SNRdB_vec,r1)
    plot(SNRdB_vec,r2)
    title( ' '.join([str(ch_index), str(sch_index), str(RX)]) )
        