# -*- coding: utf-8 -*-
"""
Created on Thu May 07 18:57:04 2015

@author: Hassan
"""

from multilayer4 import sim_init, \
generate_vecs, generate_ch_matrix, power_alloc, \
info_rate_optrx, info_rate_expand, info_rate_ild, \
save_data, load_data

#%%
from multilayer4 import SCHEMES_DICT, CHANNELS_DICT


from numpy import *
#import numpy

scheme_index = 1
channel_index = 'C'
CHANNEL_LABELS = ['A', 'B', 'C', 'D', 'E']

def run_sim(scheme_index,channel_index):
    print('Starting')
    
    PdB_vec = arange(0,30+1);
    P_vec = 10**(PdB_vec/10.0);
    
    SCH = SCHEMES_DICT[scheme_index]
    W_base, a_base, K_prime, fc_base = SCH['W_base'], SCH['a_base'], SCH['K_prime'], SCH['fc_base']
    print('Scheme parameters:')
    print('W = %f, a = %f, K^prime = %d, fc (base) = %f' % (W_base,a_base,K_prime,fc_base))
    
    CH = CHANNELS_DICT[channel_index]
    print('Channel parameters:')
    for key in ['N_paths', 'h_wb', 'tau', 'alpha']: print key + " =", CH[key]
    
    print('Initializing simulation')
    SIM_DICT = sim_init(PASSBAND=True,T_TRANSMISSION=32,F_samp=64,N0=1)
    
    print('Generating transmitter matrix')
    H_TX, f_min, f_max = generate_vecs(W_base,a_base,K_prime,fc_base,SIM_DICT)
    
    print('Generating channel matrix')
    H_CH = generate_ch_matrix( CH, SIM_DICT)
    
    print('Computing rate for optimal receiver')
    RATE_OPT = info_rate_optrx( H_TX, H_CH, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    
    print('Computing rate for expanded band receiver')
    K = K_prime + 1
    H_RX_EB, dummy, dummy = generate_vecs(W_base,a_base,K,fc_base,SIM_DICT)
    RATE_EB = info_rate_expand( H_TX, H_CH, H_RX_EB, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    
    print('Computing rate for same band receiver with joint layer decoding')
    H_RX = H_TX # the scaling is not important
    RATE_SB_JLD = info_rate_expand( H_TX, H_CH, H_RX, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    
    
    print('Computing rate for same band receiver with individual layer decoding')
    H_RX = H_TX*sqrt(2*SIM_DICT['dt']); # H_RX.transpose().dot(H_RX) = IDENTITY #XXXX                                    # the scaling is very important for correct RESULTS_PY
    RATE_SB_ILD = info_rate_ild( H_TX, H_CH, H_RX, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    
    print('Collecting results')
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
    DATA['RUNTIME'] = None
#    DATA['RATE'] = {};
#    DATA['RATE']['R_OPT'] = RATE_OPT;
#    DATA['RATE']['R_EB'] = RATE_EB;
#    DATA['RATE']['R_SB_JLD'] = RATE_SB_JLD;
#    DATA['RATE']['R_SB_ILD'] = RATE_SB_ILD;
    
    print('Finished')
    return DATA

#DATA = run_sim(scheme_index,channel_index)

#%%
def save_data_to_file(DATA):
    import datetime
    datetime_now = datetime.datetime.now()
    hms = [datetime_now.hour, datetime_now.minute, datetime_now.second]
    hms_str = ['{0:02d}'.format(i) for i in hms] # includes leading zero
    timestamp = ''.join(hms_str)
    #timestamp = str(datetime_now.hour) + str(datetime_now.minute) + str(datetime_now.second)
    ## str() ignores leading zeros
    
    tail = timestamp
    root_filename  = 'channel' + channel_index.lower() + '_' 
    root_filename += 'scheme' + str(scheme_index) + '_'
    data_filename = root_filename + tail + '_thinice_py'
    save_data(data_filename,DATA);
    print('Output saved to file: %s' %(data_filename))
    
#%%
import timeit

CHANNEL_LABELS = CHANNELS_DICT.keys()
CHANNEL_LABELS.sort()
CHANNEL_LABELS = ['A','B','E']

SCHEME_LABELS = SCHEMES_DICT.keys()
SCHEME_LABELS.sort()
SHCEME_LABELS = [1,2,3]

runtimes = list()
for channel_index in CHANNEL_LABELS:
    for scheme_index in SCHEME_LABELS:
        print('CHANNEL %s, SCHEME %d' %((channel_index,scheme_index)))
        start_time = timeit.default_timer()
    
        DATA = run_sim(scheme_index,channel_index)
            
        finish_time = timeit.default_timer()
        runtime = finish_time - start_time
        runtimes.append(runtime)
        print("Simulation runtime: %.2f sec." %(runtime))
    
        DATA['RUNTIME'] = runtime    
        save_data_to_file(DATA)
print "Done."    

#%%
#    import time
#    time.sleep(0.5) 

#%%
#RESULTS_PY = {}
#for ch_idx in range(1,5+1):
#    RESULTS_PY[('CHANNEL',ch_idx)] = {}
#    for sch_idx in range(1,3+1):
#        print("loading %s" % filename[ch_idx][sch_idx])
#        DATA = loadmat_to_dict(filename[ch_idx][sch_idx])
#        RESULTS_PY[('CHANNEL',ch_idx)][('SCHEME',sch_idx)] = DATA['SCHEME']['RATE']
#        SNRdB_vec = DATA['SNRdB'];
        