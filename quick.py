# -*- coding: utf-8 -*-
"""
Created on Wed May 13 14:32:39 2015

@author: Hassan
"""

from multilayer5 import sim_init, \
generate_vecs, generate_ch_matrix, power_alloc, \
info_rate_optrx, info_rate_expand, info_rate_ild, \
save_data, load_data, \
SCHEMES_DICT, CHANNELS_DICT

from numpy import *

#%%
def print_log(string): # prints to screen and logs to file
    print string 

#%%
scheme_index = 1
channel_index = 'C'

print_log('Starting')

PdB_vec = arange(0,30+1);
P_vec = 10**(PdB_vec/10.0);

SCH = SCHEMES_DICT[scheme_index]
W_base, a_base, K_prime, fc_base = SCH['W_base'], SCH['a_base'], SCH['K_prime'], SCH['fc_base']
print_log('Scheme parameters:')
print_log('W = %f, a = %f, K^prime = %d, fc (base) = %f' % (W_base,a_base,K_prime,fc_base))

CH = CHANNELS_DICT[channel_index]
print_log('Channel parameters:')
for key in ['N_paths', 'h_wb', 'tau', 'alpha']: print_log( key + " =" + str(CH[key]) )

print_log('Initializing simulation')
SIM_DICT = sim_init(PASSBAND=True,T_TRANSMISSION=32,F_samp=64,N0=1)
SIM_DICT = sim_init(PASSBAND=True,T_TRANSMISSION=32,F_samp=16,N0=1)

print_log('Generating transmitter matrix')
H_TX, f_min, f_max = generate_vecs(W_base,a_base,K_prime,fc_base,SIM_DICT)

print_log('Generating channel matrix')
H_CH = generate_ch_matrix( CH, SIM_DICT)

print_log('Computing rate for optimal receiver')
RATE_OPT = info_rate_optrx( H_TX, H_CH, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );

SIMRX_DICT = dict(SIM_DICT)
SIMRX_DICT['T_TRANSMISSION'] = SIMRX_DICT['T_TRANSMISSION'] + max(CH['tau'])    

print_log('Computing rate for expanded band receiver')
K = K_prime + 1
H_RX_EB, dummy, dummy = generate_vecs(W_base,a_base,K,fc_base,SIMRX_DICT)
RATE_EB = info_rate_expand( H_TX, H_CH, H_RX_EB, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );

print_log('Computing rate for same band receiver with joint layer decoding')
H_RX, dummy, dummy = generate_vecs(W_base,a_base,K_prime,fc_base,SIMRX_DICT)
RATE_SB_JLD = info_rate_expand( H_TX, H_CH, H_RX, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );

print_log('Computing rate for same band receiver with individual layer decoding')
RATE_SB_ILD = info_rate_ild( H_TX, H_CH, H_RX, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );

print_log('Collecting results')
SNRdB_vec = 10*log10(P_vec/SIM_DICT['N0'])
DATA = {};
DATA['SIM'] = SIM_DICT;
DATA['SNRdB'] = SNRdB_vec;
DATA['CHANNEL_PARAMS'] = CHANNELS_DICT[channel_index];
DATA['SCHEME_PARAMS'] = SCHEMES_DICT[scheme_index];
RX_LIST = ['OPT','EB','SB_JLD','SB_ILD']
DATA['RX'] = {}
for RX in RX_LIST: DATA['RX'][RX] = {}
DATA['RX']['OPT'   ]['RATE'] = RATE_OPT
DATA['RX']['EB'    ]['RATE'] = RATE_EB
DATA['RX']['SB_JLD']['RATE'] = RATE_SB_JLD
DATA['RX']['SB_ILD']['RATE'] = RATE_SB_ILD
DATA['RUNTIME'] = None # Simulation Runtime
DATA['RX']['EB'    ]['LABEL'] = 'EB (K=%d)' % (K) # Specify exactly how many branches were used at RX

#%%
RX_EB = {}
for K in range(1,K_prime+2):
    H_RX_EB, dummy, dummy = generate_vecs(W_base,a_base,K,fc_base,SIMRX_DICT)
    RATE_EB = info_rate_expand( H_TX, H_CH, H_RX_EB, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    RX = 'EB_K%d' % K
    RX_EB[RX] = {}
    RX_EB[RX]['RATE'] = RATE_EB
    RX_EB[RX]['LABEL'] = 'EB (K=%d)' % K
    
#%%
for RX in RX_EB.keys(): 
    DATA['RX'][RX] = {}
    DATA['RX'][RX]['RATE' ] = RX_EB[RX]['RATE' ]
    DATA['RX'][RX]['LABEL'] = RX_EB[RX]['LABEL']


#DATA = run_sim(scheme_index,channel_index)
