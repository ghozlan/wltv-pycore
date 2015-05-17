# -*- coding: utf-8 -*-
"""
Created on Thu May 07 15:27:37 2015

@author: Hassan

Reading a structure from .mat (MATLAB) file

"""

#%%

import scipy.io
from numpy import *

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

#%
def loadmat_to_dict(filename):
    mat = scipy.io.loadmat(filename)
            
    #%    
    D = parse(mat['DATA'])
    D_SNRdB = D['SNRdB']
    D_SIM = parse(D['SIM'])
    D_CHANNEL = parse(D['CHANNEL'])
    D_SCHEME = parse(D['SCHEME'])
    D_SCHEME_PARAMS = D_SCHEME['PARAMS']
    D_SCHEME_RATE = parse(D_SCHEME['RATE'])
        
    #%
    DATA ={}
    DATA['SNRdB'] = make_flat(D_SNRdB)
    DATA['SIM'] = make_flat(D_SIM)
    DATA['CHANNEL'] = make_flat(D_CHANNEL)
    DATA['SCHEME'] = {}
    DATA['SCHEME']['PARAMS'] = make_flat(D_SCHEME_PARAMS)
    DATA['SCHEME']['RATE'] = make_flat(D_SCHEME_RATE)
    
    return DATA

#%
def reformat(D): #convert matlab structure to be similar to python structure
    DATA = {};
    DATA['SIM'] = D['SIM'];
    DATA['SNRdB'] = D['SNRdB'];
    DATA['CHANNEL_PARAMS'] = D['CHANNEL']
    DATA['SCHEME_PARAMS'] = D['SCHEME']['PARAMS']
    RX_LIST = ['OPT','EB','SB_JLD','SB_ILD']
    DATA['RX'] = {}
    for RX in RX_LIST: DATA['RX'][RX] = {}
    DATA['RX']['OPT'   ]['RATE'] = D['SCHEME']['RATE']['R_OPT']
    DATA['RX']['EB'    ]['RATE'] = D['SCHEME']['RATE']['R_EB']
    DATA['RX']['SB_JLD']['RATE'] = D['SCHEME']['RATE']['R_SB_JLD']
    DATA['RX']['SB_ILD']['RATE'] = D['SCHEME']['RATE']['R_SB_ILD3'] #R_SB_ILD3 is what you want to look for (if there)
    
    return DATA
	
	
#%
def load_data(filename):
	D = loadmat_to_dict(filename)
	DATA = reformat(D)
	return DATA