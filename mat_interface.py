# -*- coding: utf-8 -*-
"""
Created on Thu May 07 15:27:37 2015

@author: Hassan

Read a structure from .mat (MATLAB) file

"""

#%%

import scipy.io
from numpy import *

#%%

# The MATLAB format
#DATA.SIM = SIM;
#DATA.SNRdB = SNRdB_vec;
#DATA.CHANNEL_PARAMS = CH;
#DATA.SCHEME_PARAMS = SCHEMES{scheme_index};
#DATA.RX.OPT   .RATE = R_OPT;
#DATA.RX.EB    .RATE = R_EB;
#DATA.RX.EB    .LABEL = sprintf('EB (K=%d)',K);
#DATA.RX.SB_JLD.RATE = R_SB_JLD;
#DATA.RX.SB_ILD.RATE = R_SB_ILD;
#DATA.RUNTIME = sim_runtime;

#%%

import os
cwd = os.getcwd()
mat_filename = cwd + '\\output-mat' + '\\' + 'channela_scheme1_003554.mat'
py_filename  = cwd + '\\output'     + '\\' + 'channela_scheme1_121140_py'

import wltvlib
pyd = wltvlib.load_data(py_filename)

#%

def parse(obj):
    keys = [obj.dtype.descr[i][0] for i in range(len(obj.dtype.descr))]
    #keys = obj.dtype.fields # this should work too
    D={}
    for i in range(len(keys)):
        D[keys[i]] = obj[0][0][i]
        #D[keys[i]] = obj[0][0][i].flatten()
    return D

#%

def make_flat(obj):
    import numpy
    if type(obj)==dict:
        new_dict = {}
        for key in obj.keys():
            new_dict[key] = obj[key].flatten()
        return new_dict
    elif type(obj)==numpy.ndarray:
        return obj.flatten()

#%%

def loadmat_to_dict(filename):
    mat = scipy.io.loadmat(filename)
            
    #%    
    D = parse(mat['DATA'])

    D_SNRdB = D['SNRdB']
    D_SIM = parse(D['SIM'])
    D_CHANNEL_PARAMS = parse(D['CHANNEL_PARAMS'])
    D_SCHEME_PARAMS = D['SCHEME_PARAMS'] #W_base, a_base, K_prime, fc_base
    
    D_RX = parse(D['RX'])

    D_RX_OPT   = parse(D_RX['OPT'   ])
    D_RX_EB    = parse(D_RX['EB'    ])
    D_RX_SBJLD = parse(D_RX['SB_JLD'])
    D_RX_SBILD = parse(D_RX['SB_ILD'])        

    D_RX_OPT_RATE   = D_RX_OPT  ['RATE']
    D_RX_EB_RATE    = D_RX_EB   ['RATE']
    D_RX_SBJLD_RATE = D_RX_SBJLD['RATE']
    D_RX_SBILD_RATE = D_RX_SBILD['RATE']    
    
    D_RX_EB_LABEL   = D_RX_EB['LABEL']
    
    #%
    DATA ={}
    DATA['SNRdB'] = make_flat(D_SNRdB)
    DATA['SIM'] = make_flat(D_SIM)
    for key in DATA['SIM'].keys(): 
        if len(DATA['SIM'][key])==1: DATA['SIM'][key] = DATA['SIM'][key][0] 
    DATA['CHANNEL_PARAMS'] = make_flat(D_CHANNEL_PARAMS)
    DATA['SCHEME_PARAMS'] = make_flat(D_SCHEME_PARAMS)

    DATA['RX'] = {}
    RX_LIST = ['OPT','EB','SB_JLD','SB_ILD']
    for RX in RX_LIST: DATA['RX'][RX] = {}
    DATA['RX']['OPT'   ]['RATE'] = make_flat(D_RX_OPT_RATE)
    DATA['RX']['EB'    ]['RATE'] = make_flat(D_RX_EB_RATE)
    DATA['RX']['EB'    ]['LABEL'] = str(make_flat(D_RX_EB_LABEL)[0])
    DATA['RX']['SB_JLD']['RATE'] = make_flat(D_RX_SBJLD_RATE)
    DATA['RX']['SB_ILD']['RATE'] = make_flat(D_RX_SBILD_RATE)
    DATA['RUNTIME'] = make_flat(D['RUNTIME'])[0]

    return DATA

#%
def reformat(D): #convert matlab structure to be similar to python structure
    W_base, a_base, K_prime, fc_base = D['SCHEME_PARAMS']
    SCH = {'W_base':W_base, 'a_base': a_base, 'K_prime': K_prime, 'fc_base':fc_base}
    
    DATA = dict(D)
    DATA['SCHEME_PARAMS'] = SCH
    
    return DATA
	
#%
def load_data(filename):
	D = loadmat_to_dict(filename)
	DATA = reformat(D)
	return DATA
 
 #%%
 
#D = load_data('channela_scheme1_003554.mat')
#D = load_data('channela_scheme2_004147.mat')
#D = load_data('channele_scheme1_002654.mat')
#D = load_data('channele_scheme2_004737.mat')