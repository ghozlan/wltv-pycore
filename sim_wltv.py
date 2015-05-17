# -*- coding: utf-8 -*-
"""
Created on Thu May 07 18:57:04 2015

@author: Hassan
"""

from multilayer5 import sim_init, \
generate_vecs, generate_ch_matrix, power_alloc, \
info_rate_optrx, info_rate_expand, info_rate_ild, \
save_data, load_data, timestamp_and_save_data, \
plot_spectrum, \
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
    sim_start = timeit.default_timer()
    
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
    SIM_DICT['T_RX'] = SIM_DICT['T_TRANSMISSION'] + max(CH['tau'])
    SNR_vec = 10**(SNRdB_vec/10.0)
    P_vec = SNR_vec * SIM_DICT['N0']
    
    print_log('Generating transmitter matrix')
    H_TX, f_min, f_max = generate_vecs(W_base,a_base,K_prime,fc_base,SIM_DICT)
    
    print_log('Generating channel matrix')
    H_CH = generate_ch_matrix( CH, SIM_DICT)

    print_log('Generating receiver matrix')
    # For same band receiver
    [H_RX, f_min, f_max] = generate_vecs(W_base,a_base,K_prime,fc_base,SIM_DICT,IS_RX=True);
    #H_RX = H_RX * SIM['dt'] # the scaling is not important
    # For expanded band receiver
    K = K_prime+1;
    [H_RX_EB, f_min, f_max] = generate_vecs(W_base,a_base,K,fc_base,SIM_DICT,IS_RX=True);
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
    start_time = timeit.default_timer()
    RATE_SB_JLD = info_rate_expand( H_TX, H_CH, H_RX, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    finish_time = timeit.default_timer()
    runtime = finish_time - start_time
    print_log('(%.2f sec)'%runtime)    
    RUNTIMES['SB_JLD'] = runtime
    
    print_log('Computing rate for same band receiver with individual layer decoding')
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
    sim_end = timeit.default_timer()
    sim_runtime = sim_end - sim_start
    DATA['RUNTIME'] = sim_runtime # Simulation Runtime
    for RX in RX_LIST: DATA['RX'][RX]['RUNTIME'] = RUNTIMES[RX]
    
    # Save data to disk    
    timestamp_and_save_data(DATA)
    
    print_log('Simulation runtime: %.2f sec' % (sim_runtime))

    # Plot input and outout spectra
    plot_spectrum(H_TX, H_CH, SCH, CH, SIM_DICT)
    fig_filename = 'SPECTRUM_CH' + channel_index +  '_SCH' + str(scheme_index)
    from matplotlib.pyplot import savefig
    savefig(fig_filename+'.eps', format='eps', dpi=1000)
    
    return DATA

#DATA = run_sim(scheme_index,channel_index)

#%% 

def run_sim_ebk(scheme_index,channel_index,SNRdB_vec=arange(0,30+1)):
    import timeit
    RUNTIMES = {}
    
    print_log('Starting')
    sim_start = timeit.default_timer()
    
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
    SIM_DICT['T_RX'] = SIM_DICT['T_TRANSMISSION'] + max(CH['tau'])
    SNR_vec = 10**(SNRdB_vec/10.0)
    P_vec = SNR_vec * SIM_DICT['N0']
    
    print_log('Generating transmitter matrix')
    H_TX, f_min, f_max = generate_vecs(W_base,a_base,K_prime,fc_base,SIM_DICT)
    
    print_log('Generating channel matrix')
    H_CH = generate_ch_matrix( CH, SIM_DICT)

    RX_EB = {}
    for K in range(1,K_prime + 2):
        EBK = 'EB_K%d' % K
    
        print_log('Generating receiver matrix')
        H_RX_EBK, dummy, dummy = generate_vecs(W_base,a_base,K,fc_base,SIM_DICT,IS_RX=True)
        #H_RX_EB = H_RX_EB * SIM['dt'] # the scaling is not important
        
        print_log('Computing rate for expanded band receiver: K = %d' %(K))
        start_time = timeit.default_timer()
        RATE_EBK = info_rate_expand( H_TX, H_CH, H_RX_EBK, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
        finish_time = timeit.default_timer()
        runtime = finish_time - start_time
        print_log('(%.2f sec)'%runtime)    
        RUNTIMES[EBK] = runtime
    
        RX_EB[EBK] = {}
        RX_EB[EBK]['RATE' ] = RATE_EBK
        RX_EB[EBK]['K'    ] = K
        RX_EB[EBK]['LABEL'] = 'EB (K=%d)' % K
    
    print_log('Collecting results')
    SNRdB_vec = 10*log10(P_vec/SIM_DICT['N0'])
    DATA = {};
    DATA['SIM'] = SIM_DICT;
    DATA['SNRdB'] = SNRdB_vec;
    DATA['CHANNEL_PARAMS'] = CHANNELS_DICT[channel_index];
    DATA['SCHEME_PARAMS'] = SCHEMES_DICT[scheme_index];
    DATA['RX'] = {}
    DATA['RX'].update(RX_EB)

    sim_end = timeit.default_timer()
    sim_runtime = sim_end - sim_start
    DATA['RUNTIME'] = sim_runtime # Simulation runtime
    RX_LIST = RX_EB.keys()
    for RX in RX_LIST: DATA['RX'][RX]['RUNTIME'] = RUNTIMES[RX]
    
    # save data to disk    
    timestamp_and_save_data(DATA,keyword='ebk')
    
    print_log('Simulation runtime: %.2f sec' % (sim_runtime))
#%%
DATA = run_sim_ebk(scheme_index,channel_index)
#%%
logfile = open('log.txt','w')
logfile.close()    

#CHANNEL_LABELS = sorted(CHANNELS_DICT.keys())
CHANNEL_LABELS = ['A','E']

#SCHEME_LABELS = sorted(SCHEMES_DICT.keys())
SCHEME_LABELS = [1, 2]

for channel_index in CHANNEL_LABELS:
    for scheme_index in SCHEME_LABELS:
        print_log('CHANNEL %s, SCHEME %d' %((channel_index,scheme_index)))        
        DATA = run_sim(scheme_index,channel_index)
        print_log( '-' * 60 )
print('Done.')