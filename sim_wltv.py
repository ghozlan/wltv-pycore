# -*- coding: utf-8 -*-
"""
Created on Thu May 07 18:57:04 2015

@author: Hassan
"""

from multilayer5 import sim_init, \
generate_vecs, generate_ch_matrix, power_alloc, \
info_rate_optrx, info_rate_expand, info_rate_ild, \
save_data, load_data, timestamp_and_save_data, \
SCHEMES_DICT, CHANNELS_DICT

from numpy import *

#%%
def print_log(string): # prints to screen and logs to file
    print string 
    with open('log.txt','a') as logfile:
        logfile.write(string + '\n')

#%%
scheme_index = 1
channel_index = 'C'

def run_sim(scheme_index,channel_index,SNRdB_vec=arange(0,30+1)):
    import timeit
    RUNTIMES = {}
    
    print_log('Starting')
    
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
    SNR_vec = 10**(SNRdB_vec/10.0)
    P_vec = SNR_vec * SIM_DICT['N0']
    
    print_log('Generating transmitter matrix')
    H_TX, f_min, f_max = generate_vecs(W_base,a_base,K_prime,fc_base,SIM_DICT)
    
    print_log('Generating channel matrix')
    H_CH = generate_ch_matrix( CH, SIM_DICT)

    print_log('Generating receiver matrix')
    # Make receiver take into account the max delay (seems to be not critical)
    SIM_RX = dict(SIM_DICT);
    SIM_RX['T_TRANSMISSION'] = SIM_DICT['T_TRANSMISSION'] + max(CH['tau']);
    # For same band receiver
    [H_RX, f_min, f_max] = generate_vecs(W_base,a_base,K_prime,fc_base, SIM_RX);
    # For expanded band receiver
    K = K_prime+1;
    [H_RX_EB, f_min, f_max] = generate_vecs(W_base,a_base,K,fc_base, SIM_RX);
    #H_RX_EB = H_RX_EB * SIM['dt'] # the scaling is not important

    
    print_log('Computing rate for optimal receiver')
    start_time = timeit.default_timer()        
    RATE_OPT = info_rate_optrx( H_TX, H_CH, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    finish_time = timeit.default_timer()
    runtime = finish_time - start_time
    print_log('(%.2f sec)'%runtime)
    RUNTIMES['OPT'] = runtime
    
    print_log('Computing rate for expanded band receiver')
    start_time = timeit.default_timer()
    RATE_EB = info_rate_expand( H_TX, H_CH, H_RX_EB, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    finish_time = timeit.default_timer()
    runtime = finish_time - start_time
    print_log('(%.2f sec)'%runtime)    
    RUNTIMES['EB'] = runtime
    
    print_log('Computing rate for same band receiver with joint layer decoding')
    #H_RX = H_TX # the scaling is not important
    H_RX, dummy, dummy = generate_vecs(W_base,a_base,K_prime,fc_base,SIMRX_DICT)
    start_time = timeit.default_timer()
    RATE_SB_JLD = info_rate_expand( H_TX, H_CH, H_RX, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    finish_time = timeit.default_timer()
    runtime = finish_time - start_time
    print_log('(%.2f sec)'%runtime)    
    RUNTIMES['SB_JLD'] = runtime
    
    print_log('Computing rate for same band receiver with individual layer decoding')
    #H_RX = H_TX*sqrt(2*SIM_DICT['dt']); # H_RX.transpose().dot(H_RX) = IDENTITY #XXXX                                    
    # the scaling is very important for correct RESULTS_PY
    start_time = timeit.default_timer()
    RATE_SB_ILD = info_rate_ild( H_TX, H_CH, H_RX, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    finish_time = timeit.default_timer()
    runtime = finish_time - start_time
    print_log('(%.2f sec)'%runtime)    
    RUNTIMES['SB_ILD'] = runtime
    
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
    DATA['RX']['EB'    ]['LABEL'] = 'EB (K=%d)' % (K)
    DATA['RX']['SB_JLD']['RATE'] = RATE_SB_JLD
    DATA['RX']['SB_ILD']['RATE'] = RATE_SB_ILD
    DATA['RUNTIME'] = None # Simulation Runtime
    for RX in RX_LIST: DATA['RX'][RX]['RUNTIME'] = RUNTIMES[RX]
    
    RX_EB = {}
    #for K in range(1,K_prime + 2):
    H_RX_EBK, dummy, dummy = generate_vecs(W_base,a_base,K,fc_base,SIMRX_DICT)
    RATE_EBK = RATE_EB
    #RATE_EBK = info_rate_expand( H_TX, H_CH, H_RX_EB, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    EBK = 'EB_K%d' % K
    RX_EB[EBK] = {}
    RX_EB[EBK]['RATE' ] = RATE_EBK
    RX_EB[EBK]['LABEL'] = 'EB (K=%d)' % K
    
    print('Finished')
    return DATA

#DATA = run_sim(scheme_index,channel_index)

#%%
def timestamp_and_save_data(DATA,keyword=''):
    import datetime
    datetime_now = datetime.datetime.now()
    hms = [datetime_now.hour, datetime_now.minute, datetime_now.second]
    hms_str = ['{0:02d}'.format(i) for i in hms] # includes leading zero
    timestamp = ''.join(hms_str)
    #timestamp = str(datetime_now.hour) + str(datetime_now.minute) + str(datetime_now.second)
    ## str() ignores leading zeros
    
    tail = timestamp + '_py' + keyword
    root_filename  = 'channel' + channel_index.lower() + '_' 
    root_filename += 'scheme' + str(scheme_index) + '_'
    data_filename = root_filename + tail
    save_data(data_filename,DATA);
    print('Output saved to file: %s' %(data_filename))

#timestamp_and_save_data(DATA)    
#%%
logfile = open('log.txt','w')
logfile.close()    

CHANNEL_LABELS = CHANNELS_DICT.keys()
CHANNEL_LABELS.sort()
CHANNEL_LABELS = ['A','E']

SCHEME_LABELS = SCHEMES_DICT.keys()
SCHEME_LABELS.sort()
SCHEME_LABELS = [1, 2]

import timeit
runtimes = list()
for channel_index in CHANNEL_LABELS:
    for scheme_index in SCHEME_LABELS:
        print_log('CHANNEL %s, SCHEME %d' %((channel_index,scheme_index)))
        start_time = timeit.default_timer()    
        
        DATA = run_sim(scheme_index,channel_index)            

        finish_time = timeit.default_timer()
        runtime = finish_time - start_time
        DATA['RUNTIME'] = runtime    
        print_log('Simulation runtime: %.2f sec.' %(runtime))
        print_log( '-' * 60 )
        runtimes.append(runtime)
        #timestamp_and_save_data(DATA)
print('Done.')