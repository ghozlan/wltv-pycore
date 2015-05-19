# -*- coding: utf-8 -*-
"""
Created on Thu May 07 18:57:04 2015

Simulations

@author: Hassan
"""

from wltvlib import SCHEMES_DICT, CHANNELS_DICT
from multilayer import *

from numpy import *

#%% Set signal-to-noise ratio (SNR)
SNRdB_vec = arange(0,30+1);

#%% Channel A, Scheme 1

channel_index = 'A'
scheme_index = 1
print('CHANNEL %s, SCHEME %d' %((channel_index,scheme_index)))        
run_sim(scheme_index,channel_index,SNRdB_vec)
print( '-' * 60 )

#%% Channel A, Scheme 2

channel_index = 'A'
scheme_index = 2
print('CHANNEL %s, SCHEME %d' %((channel_index,scheme_index)))        
run_sim(scheme_index,channel_index,SNRdB_vec)
print( '-' * 60 )

#%% Channel E, Scheme 2

channel_index = 'E'
scheme_index = 2
print('CHANNEL %s, SCHEME %d' %((channel_index,scheme_index)))        
run_sim(scheme_index,channel_index,SNRdB_vec)
print( '-' * 60 )


##%% Loop
#
## Set signal-to-noise ratio (SNR)
#SNRdB_vec = arange(0,30+1);
#
##CHANNEL_LABELS = sorted(CHANNELS_DICT.keys())
#CHANNEL_LABELS = ['A','E']
#
##SCHEME_LABELS = sorted(SCHEMES_DICT.keys())
#SCHEME_LABELS = [1, 2]
#
#for channel_index in CHANNEL_LABELS:
#    for scheme_index in SCHEME_LABELS:
#        print('CHANNEL %s, SCHEME %d' %((channel_index,scheme_index)))        
#        DATA = run_sim(scheme_index,channel_index,SNRdB_vec)
#        print_log( '-' * 60 )
#print('Done.')