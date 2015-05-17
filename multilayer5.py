# -*- coding: utf-8 -*-
"""
Created on Tue May 05 17:27:14 2015
# USE DICTIONARY RATHER THAN CLASS
# LOOP
@author: Hassan
"""


from __future__ import division #makes float division default. for integer division, use //. e.g. 7//4=1

from numpy import *
from matplotlib.pyplot import *
from numpy.fft import ifft, fftshift, fft
from numpy.linalg import svd
from numpy.linalg import eigvalsh #eigvalsh : eigenvalues of a symmetric or Hermitian (conjugate symmetric)

#%%

if __name__ == '__main__':
    print "Running as main"
else:
    print "Importing module"

#%%
def print_log(string): # prints to screen and logs to file
    print string 
    
#%% Choose Scheme, Channel, and SNR values

scheme_index = 1        # 1, 2, 3, 4
channel_index = 'E';    # 'A', 'B', 'C', 'D', 'E'

SNRdB_vec = arange(0,30+1);
RATE_ZEROS = zeros(SNRdB_vec.shape)

#%% Simulation Parameters
def sim_init(PASSBAND=True,T_TRANSMISSION=32,F_samp=64,N0=1):
    SIM = {}
    SIM['PASSBAND'] = PASSBAND;
    if(PASSBAND):
        SIM['REAL_DIM_PER_SYM'] = 1;
    else:
        SIM['REAL_DIM_PER_SYM'] = 2;

    SIM['T_TRANSMISSION'] = T_TRANSMISSION;
    SIM['T_SIMULATION'] = 3*T_TRANSMISSION;        
    SIM['df'] = 1.0/SIM['T_SIMULATION'];
    
    SIM['F_samp'] = F_samp;
    SIM['dt'] = 1.0/F_samp;

    t = linspace(0,SIM['T_SIMULATION'],SIM['T_SIMULATION']/SIM['dt']+1); #t = 0:dt:T_SIMULATION; 
    t = t[:-1] # t=t(1:end-1), which yields same effect as: t(end) = [];
    #t = delete(t, len(t)-1)     #t(end) = [];
    SIM['t'] = t
    
    SIM['N0'] = N0;

    return SIM


SIM_DICT = sim_init(PASSBAND=True,T_TRANSMISSION=32,F_samp=16,N0=1)

PASSBAND = SIM_DICT['PASSBAND'];
REAL_DIM_PER_SYM = SIM_DICT['REAL_DIM_PER_SYM'];
T_TRANSMISSION = SIM_DICT['T_TRANSMISSION'];
T_SIMULATION = SIM_DICT['T_SIMULATION'];
df = SIM_DICT['df'];
F_samp = SIM_DICT['F_samp'];
dt = SIM_DICT['dt'];
t = SIM_DICT['t'];
N0 = SIM_DICT['N0'];

SNR_vec = 10**(SNRdB_vec/10.0);
P_vec = SNR_vec * SIM_DICT['N0']

#%% SCHEMES

SCHEMES_DICT = {}
SCHEMES_DICT[1] = {'W_base':1, 'a_base': 2, 'K_prime': 3, 'fc_base':1.5}
SCHEMES_DICT[2] = {'W_base':1, 'a_base': 1.587401051968199, 'K_prime': 4, 'fc_base':1.5}
SCHEMES_DICT[3] = {'W_base':7, 'a_base': 1, 'K_prime': 1, 'fc_base':4.5}
SCHEMES_DICT[4] = {'W_base':1, 'a_base': 1, 'K_prime': 1, 'fc_base':1.5}

if __name__ == '__main__':
    SCH = SCHEMES_DICT[scheme_index]
    W_base, a_base, K_prime = SCH['W_base'], SCH['a_base'], SCH['K_prime']
    fc_base = SCH['fc_base']
    
    print('Scheme parameters:')
    print('W = %f, a = %f, K^prime = %d, fc (base) = %f' % (W_base,a_base,K_prime,fc_base))

#%% Generating Transmitter Matrix H_TX
def generate_vecs(W_base,a_base,N_layers,fc_base, SIM,rx=False):
    # return [H_TX, f_min, f_max]

    PASSBAND = SIM['PASSBAND']
    T_XX = SIM['T_TRANSMISSION']
    if rx==True:
        T_XX = SIM['T_RX']    
    F_samp = SIM['F_samp']
    T_SIMULATION = SIM['T_SIMULATION']
    dt = SIM['dt']
    df = SIM['df']

    
    a_vec = a_base**arange(0,N_layers);
    
    T_sym = 1/W_base;
    N_sym = int(T_XX/T_sym);  # number of signals (in base layer) #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    L = int(T_sym/dt);             # oversampling factor (for base layer)


    H = {};
    k_prime = 1;
    for a in a_vec:
        W = W_base * a;         # bandwidth
        
        if(PASSBAND):
            fc = fc_base * a;       # carrier frequency
        else:
            fc = 0;
        
     
        N_FFT = T_SIMULATION/dt;
        f = linspace(-F_samp/2,F_samp/2,N_FFT);
        Xf = zeros(len(f));
        Xf[abs(f-fc)<W/2] = 1/sqrt(W);
        Xf = fftshift(Xf);
        xt = fftshift(N_FFT*ifft(Xf*df));
        xt_pb_iphase =  sqrt(2)*real(xt); #pass band: sum(abs(xt_pb_iphase).^2)*dt = 1
        xt_pb_qphase = -sqrt(2)*imag(xt); #pass band: sum(abs(xt_pb_qphase).^2)*dt = 1
        
        
        N = int(floor(N_sym * a));
        #H_TX = zeros(length(xt),N); % tall matrix
        
        if(PASSBAND):
            H_TX_IPHASE = zeros((len(xt),N)); # tall matrix
            H_TX_QPHASE = zeros((len(xt),N)); # tall matrix
            for n in range(N):
                H_TX_IPHASE[:,n] = roll(xt_pb_iphase,int( (n-1)*ceil(L/a) )); # passband inphase
                H_TX_QPHASE[:,n] = roll(xt_pb_qphase,int( (n-1)*ceil(L/a) )); # passband quadrature phase
                
            H_TX = 1/sqrt(2) * hstack((H_TX_IPHASE,H_TX_QPHASE)); # inphase energy = 1/2, quadrature phase energy = 1/2
            print('Layer %d: 2x%d (real) symbols' %(k_prime,N))
        else:
            H_TX_BASEBAND = zeros((len(xt),N)); # tall matrix
            for n in range(N):
                H_TX_BASEBAND[:,n] = roll(xt,int( (n-1)*ceil(L/a) ));
            
            H_TX = H_TX_BASEBAND;
            print('Layer %d: %d (complex) symbols' %(k_prime,N))
        
        
        H[k_prime] = H_TX; # H{i} is a matrix in which the columns are the tx vector of layer i
        k_prime = k_prime + 1;
    
    f_min = min( a_vec*(fc_base-W_base/2) );
    f_max = max( a_vec*(fc_base+W_base/2) );
    B = f_max-f_min;
    print('f_min = %.5f, f_max = %.5f, B = %.5f' %(f_min,f_max,B))
    
    
    #K_prime = N_layers;
    #%B = W_base * (a_base.^K_prime-1)/(a_base-1) % THIS EXPRESSION IS RIGHT
    #ONLY WHEN LAYERS ARE ORTHOGONAL IN FREQUENCY
    
    #Concatenate all tx (column) vectors in one matrix
    H_TX = H[1];
    for k_prime in range(2,N_layers+1):
        H_TX = hstack( (H_TX, H[k_prime]) );

    return [H_TX, f_min, f_max]
    #return H_TX
    
#%% 

if __name__ == '__main__':
    print_log('Generating transmitter matrix')
    H_TX, f_min, f_max = generate_vecs(W_base,a_base,K_prime,fc_base,SIM_DICT)
    B_TOTAL = f_max - f_min
    f_center = (f_max + f_min) / 2.0

#%% CHANNELS

CHANNELS_DICT = {}
CHANNELS_DICT['A'] = {'N_paths':2, 'h_wb':[1, 1/2], 'tau':[0, 2], 'alpha':[1, 2]}
CHANNELS_DICT['B'] = {'N_paths':2, 'h_wb':[1, 1/2], 'tau':[0, 2], 'alpha':[1, 1.587401051968199]}
CHANNELS_DICT['C'] = {'N_paths':2, 'h_wb':[1, 1.5], 'tau':[0, 2], 'alpha':[1, 2]}
CHANNELS_DICT['D'] = {'N_paths':2, 'h_wb':[1, 1.5], 'tau':[2, 3], 'alpha':[1, 2]}
CHANNELS_DICT['E'] = {'N_paths':3, 'h_wb':[1, -0.7, 1.5],'tau':[2, 1, 3],'alpha':[1, 1.25, 2]}

if __name__ == '__main__':
    CH = CHANNELS_DICT[channel_index]
    
    print('Channel parameters:')
    for key in ['N_paths', 'h_wb', 'tau', 'alpha']: print key + " =", CH[key]
    
    SIM_DICT['T_RX'] = SIM_DICT['T_TRANSMISSION'] + max(CH['tau']);    
        
#%% Generating the channel matrix H_CH

def generate_ch_matrix( CH, SIM ):
    # return K0_t_tau
    # Supports Widebandband Linear Time-Varying (LTV) Channels

    t, dt, F_samp = SIM['t'], SIM['dt'], SIM['F_samp'];

    N_paths, h_wb, tau, alpha = CH['N_paths'], CH['h_wb'], CH['tau'], CH['alpha']
    
    [tau_,t_] = meshgrid(t,t);
    h_t_tminustau = zeros(len(t));
    for m in range(N_paths):
        h_t_tminustau = h_t_tminustau + h_wb[m] * sqrt(alpha[m]) * \
         (1/dt) * sinc(F_samp*( (t_-tau_) - (alpha[m]*tau[m]-(alpha[m]-1)*t_) ));
    
    K0_t_tau = h_t_tminustau; # Kernel
    
    return K0_t_tau * dt

#%%

if __name__ == '__main__':
    print_log('Generating channel matrix')
    H_CH = generate_ch_matrix( CH, SIM_DICT)

#%% 

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
    
    
#%% Generate covariance matrix of the transmitted symbols
#   and some information about the layers 

def power_alloc( H_TX, SCHEME, SIM ):
#return [Sigma_X_NORMALIZED, layer]

    dt = SIM['dt']; T_TRANSMISSION = SIM['T_TRANSMISSION']; PASSBAND = SIM['PASSBAND'];
    
    W_base = SCHEME['W_base']; a_base = SCHEME['a_base']; K_prime = SCHEME['K_prime']    
    
    
    W_vec = W_base * a_base**arange(K_prime)
    N_symb_per_layer = floor(T_TRANSMISSION * W_vec) #XXXXXXXXXXXXXXXXXXXXXXXXXXXX
    if(PASSBAND): N_symb_per_layer *= 2
    N_sym_total = sum(N_symb_per_layer)
    
    #print N_symb_per_layer
    #N_symb_per_layer = N_symb_per_layer.astype(numpy.int64) # integer
    N_symb_per_layer = N_symb_per_layer.astype(int64) # integer
    #print N_symb_per_layer
    
    #first_index_of_layer = cumsum([0 N_symb_per_layer(1:end-1)]) + 1
    first_index_of_layer = cumsum( hstack((0,N_symb_per_layer[0:-1])) ) + 1
    last_index_of_layer = cumsum( N_symb_per_layer[0:] )
    
    first_index_of_layer -= 1
    last_index_of_layer -= 1
    
    #print first_index_of_layer
    #print last_index_of_layer
    
    Sigma_X_DIAG = ones(N_sym_total); 
    LAYER_MASK = zeros((N_sym_total,N_sym_total));
    layer = {};
    for layer_index in range(K_prime):
        layer[layer_index] = range(first_index_of_layer[layer_index],last_index_of_layer[layer_index]+1);
        Sigma_X_DIAG[layer[layer_index]] = 1/W_vec[layer_index];
        #LAYER_MASK[layer[layer_index],layer[layer_index]] = 1;
        for j in layer[layer_index]: LAYER_MASK[layer[layer_index],j] = 1;
    
    Sigma_X_NORMALIZED = 1/K_prime * diag(Sigma_X_DIAG);
    
    POWER_FACTOR = trace((H_TX.transpose().dot(H_TX).dot(Sigma_X_NORMALIZED)))*dt/T_TRANSMISSION #SHOULD BE = 1
    Sigma_X_NORMALIZED = Sigma_X_NORMALIZED / POWER_FACTOR;
    
    #print('CHECK POWER CONSTRAINT: POWER = %f (linear scale)' % POWER_FACTOR)
    
    return [Sigma_X_NORMALIZED, layer]

#%%

if __name__ == '__main__':
    #Sigma_X_NORMALIZED, layer = power_alloc( H_TX, SCHEMES_DICT[scheme_index], SIM_DICT)
    pass
    
#%% Rate of optimal receiver
    
def info_rate_optrx( H_TX, H_CH, P_vec, SCHEME, SIM ):
    #returns R_vec

    N0 = SIM['N0'];
    dt = SIM['dt'];
    REAL_DIM_PER_SYM = SIM['REAL_DIM_PER_SYM'];
    T_TRANSMISSION = SIM['T_TRANSMISSION'];
    
    #%--------------------------------------------------------
    N_sym_total = int(H_TX.shape[1]) #N_sym_total = size(H_TX,2);
    
    Sigma_X_NORMALIZED, layer = power_alloc(H_TX, SCHEME, SIM);
    #Sigma_X_NORMALIZED =  eye(N_sym_total); #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    POWER_FACTOR = trace((H_TX.transpose().dot(H_TX).dot(Sigma_X_NORMALIZED)))*dt/T_TRANSMISSION #SHOULD BE = 1
    Sigma_X_NORMALIZED = Sigma_X_NORMALIZED / POWER_FACTOR;
    
    #%--------------------------------------------------------
    H_CHTX = H_CH.dot(H_TX);
    #H_CHTX = H_CH.dot(H_TX);
    N_sim = int(H_TX.shape[0])       #N_sim = size(H_TX,1);
    
    sigma2_N = N0/2 * REAL_DIM_PER_SYM * (1/dt);  # DEPENDS ON WHETHER NOISE IS REAL (PASSBAND) OR COMPLEX (BASEBAND)
    #trace(Sigma_N)*dt
    
    # Sigma_Z = Sigma_N * eye(NN)
    HYX = N_sim * log(sigma2_N);
    
    Sigma_S_NORMALIZED = H_CHTX.dot(Sigma_X_NORMALIZED).dot(H_CHTX.transpose()); 
    EIG_Sigma_S_NORMALIZED = eigvalsh(Sigma_S_NORMALIZED);
    
    EIG_Sigma_Z = sigma2_N * ones(N_sim);
    #%--------------------------------------------------------
    I_vec = zeros(len(P_vec));
    loop_index = 0;
    for P in P_vec:
        #Sigma_X = P * Sigma_X_NORMALIZED;
        
        EIG_Sigma_Y = P * EIG_Sigma_S_NORMALIZED + EIG_Sigma_Z; 
        
        HY = sum(log( EIG_Sigma_Y ));
        I = HY - HYX;
        I_vec[loop_index] = I;
        loop_index = loop_index + 1;
    
    
    R_vec = I_vec / 2 * REAL_DIM_PER_SYM / T_TRANSMISSION;
    
    return R_vec

#%%

if __name__ == '__main__':
    print_log('Computing rate for optimal receiver')
    RATE_OPT = info_rate_optrx( H_TX, H_CH, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );

#%% Rate of expanded band receiver with joint decoding

def info_rate_expand( H_TX, H_CH, H_RX, P_vec, SCHEME, SIM ):
    #return R_vec

    #[U, S, V] = svd(H_RX);
    #svalue = diag(S);
    #SELECTED_COLS_IDX = find(svalue > max(svalue) * 1e-4);
    #H_RX_REDUCED = U(:,SELECTED_COLS_IDX);
    
    [U, svalue, V] = svd(H_RX);
    IS_ABOVE_THRESHOLD = svalue > max(svalue) * 1e-4;
    SELECTED_COLS_IDX = [i for i in range(len(IS_ABOVE_THRESHOLD)) if IS_ABOVE_THRESHOLD[i]==True];
    H_RX_REDUCED = U[:,SELECTED_COLS_IDX];
    
    #%-------------------------------------------------------
    N0 = SIM['N0'];
    dt = SIM['dt'];
    REAL_DIM_PER_SYM = SIM['REAL_DIM_PER_SYM'];
    T_TRANSMISSION = SIM['T_TRANSMISSION'];
    
    
    Sigma_X_NORMALIZED, layer = power_alloc(H_TX, SCHEME, SIM);
    
    H_RXCHTX = H_RX_REDUCED.transpose().dot(H_CH).dot(H_TX);
    
    sigma2_N = N0/2 * REAL_DIM_PER_SYM * (1/dt);  # DEPENDS ON WHETHER NOISE IS REAL (PASSBAND) OR COMPLEX (BASEBAND)
    #trace(Sigma_N)*dt
    
    Sigma_Z = H_RX_REDUCED.transpose().dot(sigma2_N).dot(H_RX_REDUCED);
    trace(Sigma_Z)*dt
    HYX = sum(log(eigvalsh(Sigma_Z))) #H(Y|X)
        
    I_vec = zeros(len(P_vec));
    loop_index = 0;
    for P in P_vec:
        Sigma_X = P * Sigma_X_NORMALIZED; #XXXXXXXXXXXXXXXXXXXXXXXXXXXX
        Sigma_S = H_RXCHTX.dot(Sigma_X).dot(H_RXCHTX.transpose()); 
        HY = sum(log(eigvalsh(Sigma_S + Sigma_Z)))
        I = HY - HYX;
        I_vec[loop_index] = I;
        loop_index = loop_index + 1;
    
    R_vec = I_vec / 2 * REAL_DIM_PER_SYM / T_TRANSMISSION;

    return R_vec
#%%    

if __name__ == '__main__':
    print_log('Computing rate for expanded band receiver')
    RATE_EB = info_rate_expand( H_TX, H_CH, H_RX_EB, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );

#%% Rate of same band receiver with joint layer decoding

#H_RX = H_TX * SIM['dt'];
if __name__ == '__main__':
    print_log('Computing rate for same band receiver with joint layer decoding')
    H_RX = H_TX # the scaling is not important
    RATE_SB_JLD = info_rate_expand( H_TX, H_CH, H_RX, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );

#%% Rate of same band receiver with individual layer decoding

def info_rate_ild( H_TX, H_CH, H_RX, P_vec, SCHEME, SIM ):
    N0 = SIM['N0'];
    dt = SIM['dt'];
    REAL_DIM_PER_SYM = SIM['REAL_DIM_PER_SYM'];
    T_TRANSMISSION = SIM['T_TRANSMISSION'];
    
    # Power allocation (uniform)
    Sigma_X_NORMALIZED, layer = power_alloc(H_TX, SCHEME, SIM);
    K_prime = len(layer); #K_prime = SHCEME['K_prime']
    
    # Compute mutual informal layer by layer
    R_vec_per_layer = zeros((K_prime,len(P_vec)));
    for k in range(K_prime):
        sigma2_N = N0/2 * REAL_DIM_PER_SYM * (1/dt);  # DEPENDS ON WHETHER NOISE IS REAL (PASSBAND) OR COMPLEX (BASEBAND)
        #Sigma_N = sigma2_N * eye(length(layer{k}));
        
        H_RX_k = H_RX[:,layer[k]];
        H_RXCHTX = H_RX_k.transpose().dot(H_CH).dot(H_TX);

        Sigma_Z = H_RX_k.transpose().dot(sigma2_N).dot(H_RX_k);        
            
        #EIG_S = eigvalsh(Sigma_S);
        #EIG_S_bar = eigvalsh(Sigma_S_bar);
        
        I_vec = zeros(len(P_vec));
        loop_index = 0;
        for P in P_vec:
            Sigma_X = P * Sigma_X_NORMALIZED;     #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
            Sigma_S = H_RXCHTX.dot(Sigma_X).dot(H_RXCHTX.transpose());
            
            Sigma_X_bar = Sigma_X;
            Sigma_X_bar[:,layer[k]] = 0;
            Sigma_X_bar[layer[k],:] = 0;
            
            Sigma_S_bar = H_RXCHTX.dot(Sigma_X_bar).dot(H_RXCHTX.transpose());

            I = sum(log(eigvalsh( Sigma_S+Sigma_Z ))) - sum(log(eigvalsh( Sigma_S_bar+Sigma_Z )));
            #I = sum(log( P * EIG_S + sigma2_N )) - sum(log( P * EIG_S_bar + sigma2_N ));
            I_vec[loop_index] = I;
            loop_index = loop_index + 1;
        
        R_vec_per_layer[k,:] = I_vec / 2 * REAL_DIM_PER_SYM / T_TRANSMISSION;

    
    R_vec = sum(R_vec_per_layer,0); # sum the rows (sum along the columns)
    
    return R_vec
    
#%% Rate of same band receiver with individual layer decoding ("fast" implementation) 
# turns out to be slower

def info_rate_ild_fast( H_TX, H_CH, H_RX, P_vec, SCHEME, SIM ):
    N0 = SIM['N0'];
    dt = SIM['dt'];
    REAL_DIM_PER_SYM = SIM['REAL_DIM_PER_SYM'];
    T_TRANSMISSION = SIM['T_TRANSMISSION'];
    
    # Power allocation (uniform)
    Sigma_X_NORMALIZED, layer = power_alloc(H_TX, SCHEME, SIM);
    K_prime = len(layer); #K_prime = SHCEME['K_prime']
    
    # Compute mutual informal layer by layer
    R_vec_per_layer = zeros((K_prime,len(P_vec)));
    for k in range(K_prime):
        sigma2_N = N0/2 * REAL_DIM_PER_SYM * (1/dt);  # DEPENDS ON WHETHER NOISE IS REAL (PASSBAND) OR COMPLEX (BASEBAND)
        #Sigma_N = sigma2_N * eye(length(layer{k}));
        
        H_RX_k = H_RX[:,layer[k]];
        H_RX_k, s, V = svd(H_RX_k); # ensure that H_RX_k is orthonormal
        H_RXCHTX = H_RX_k.transpose().dot(H_CH).dot(H_TX);
        
        
        Sigma_X = Sigma_X_NORMALIZED;     #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        Sigma_S = H_RXCHTX.dot(Sigma_X).dot(H_RXCHTX.transpose());
        
        #Sigma_Z = H_RX_k.transpose().dot(sigma2_N).dot(H_RX_k);
        
        Sigma_X_bar = Sigma_X;
        Sigma_X_bar[:,layer[k]] = 0;
        Sigma_X_bar[layer[k],:] = 0;
        
        Sigma_S_bar = H_RXCHTX.dot(Sigma_X_bar).dot(H_RXCHTX.transpose());
    
        EIG_S = eigvalsh(Sigma_S);
        EIG_S_bar = eigvalsh(Sigma_S_bar);
        
        I_vec = zeros(len(P_vec));
        loop_index = 0;
        for P in P_vec:
            #I = sum(log(eig( P*Sigma_S+Sigma_Z ))) - sum(log(eig( P*Sigma_S_bar+Sigma_Z )));
            I = sum(log( P * EIG_S + sigma2_N )) - sum(log( P * EIG_S_bar + sigma2_N ));
            I_vec[loop_index] = I;
            loop_index = loop_index + 1;
        
        R_vec_per_layer[k,:] = I_vec / 2 * REAL_DIM_PER_SYM / T_TRANSMISSION;

    
    R_vec = sum(R_vec_per_layer,0); # sum the rows (sum along the columns)
    
    return R_vec
    
#%%

if __name__ == '__main__':
    print_log('Computing rate for same band receiver with individual layer decoding')

    H_RX = H_TX*sqrt(2*SIM_DICT['dt']); # H_RX.transpose().dot(H_RX) = IDENTITY #XXXX
                                        # the scaling is very important for correct results
    RATE_SB_ILD = info_rate_ild( H_TX, H_CH, H_RX, P_vec, SCHEMES_DICT[scheme_index], SIM_DICT );
    

#%%
import sys
#if __name__ == '__main__': sys.exit("Finished.")
    
#%% Collect the results of the simulation
    
if __name__ == '__main__': 
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
    
    

#%% Save the results to file

def save_data(filename,py_object):
    import pickle
    with open(filename, 'w') as fout: #no need to close the file
      pickle.dump(py_object, fout)

def load_data(filename):
    import pickle
    with open(filename, 'r') as fin: #no need to close the file
      py_object = pickle.load(fin)
    return py_object
    
def timestamp_and_save_data(DATA,keyword=''):
    import datetime
    datetime_now = datetime.datetime.now()
    hms = [datetime_now.hour, datetime_now.minute, datetime_now.second]
    hms_str = ['{0:02d}'.format(i) for i in hms] # includes leading zero
    timestamp = ''.join(hms_str)
    #timestamp = str(datetime_now.hour) + str(datetime_now.minute) + str(datetime_now.second)
    ## str() ignores leading zeros
    
    tail = timestamp + '_' + keyword + '_py'
    root_filename  = 'channel' + channel_index.lower() + '_' 
    root_filename += 'scheme' + str(scheme_index) + '_'
    data_filename = root_filename + tail
    save_data(data_filename,DATA);
    print('Output saved to file: %s' %(data_filename))
    

if __name__ == '__main__': 
    timestamp_and_save_data(DATA,keyword='wahran')
    
#%%

if __name__ == '__main__':     
    #save_data('multilayer.pickle',DATA_DICT)
    #D_DICT = load_data('multilayer.pickle')
    pass


#%% Plot Input and Outout Spectra (broken down by layer)

def plot_spectrum(H_TX, H_CH, SCH, CH, SIM):
    from numpy.fft import ifft, fftshift, fft
    [Sigma_X_NORMALIZED, layer] = power_alloc(H_TX, SCH, SIM_DICT);
    K_prime = len(layer);
    SELECTED_VECS = list()
    for k in range(K_prime):
        SELECTED_VECS.append(layer[k][0]);
        
        
    N_FFT = SIM['T_SIMULATION']/SIM['dt'];
    f = linspace(-SIM['F_samp']/2,SIM['F_samp']/2,N_FFT);
    
    H_TX_SELECTED = H_TX[:,SELECTED_VECS];
    FT_IN = fft(H_TX_SELECTED, axis=0)*SIM['dt'];
    FT_IN = fftshift(FT_IN, axes=0)
    
    H_CHTX = H_CH.dot(H_TX_SELECTED);
    FT_OUT = fft(H_CHTX, axis=0)*SIM['dt'];
    FT_OUT = fftshift(FT_OUT, axes=0);
    
    #
    W_base, a_base, K_prime, fc_base = SCH['W_base'], SCH['a_base'], SCH['K_prime'], SCH['fc_base']
    f_max = (fc_base+W_base/2.0)*(a_base**(K_prime-1));
    F = f_max * max(CH['alpha']);
    
    #%--------------------------------------------------
    LineWidth = 2;
    MarkerSize = 10;
    FontSize = 17;
    FontSizeLegend = 14;
    
    line_color = ['b',[0, 0.5, 0],'r','m'];
    
    from matplotlib.pyplot import figure, subplot, plot, xlabel, ylabel, title, show
    
    fig = figure()
    subplot(K_prime+1,1,1)
    for k in range(K_prime):
        line, = plot(f,abs(FT_IN[:,k]), color=line_color[k], linewidth=LineWidth)
        #axis([-F F 0 1])
    xlabel('f')
    ylabel('|FT|')
    title('Input Spectrum')
    
    for k in range(K_prime):
        subplot(K_prime+1,1,k+2)
        plot(f,abs(FT_OUT[:,k]), color=line_color[k], linewidth=LineWidth)
        #axis([-F F 0 1])
        xlabel('f')
        ylabel('|FT|')
        title( 'Output Spectrum (k''=%d)'%k )
    
    return fig

#set(20,'Position',[2*480 50 2*480 2*470]);figure(20)

if __name__ == '__main__':
    plot_spectrum(H_TX, H_CH, SCH, CH, SIM_DICT)
    fig_filename = 'SPECTRUM_CH' + channel_index +  '_SCH' + str(scheme_index)
    savefig(fig_filename+'.eps', format='eps', dpi=1000)
    savefig(fig_filename+'.png')
    import Image
    Image.open(fig_filename+'.png').save(fig_filename+'.jpg','JPEG')


#%%
#SNRdB_vec = 10*log10(P_vec/N0)
#figure(10)
#plot(SNRdB_vec,RATE_OPT/log(2), 'b-o', label='Opt')
#plot(SNRdB_vec,RATE_SB_JLD/log(2), 'r-x', label='SB JLD')
#plot(SNRdB_vec,RATE_SB_ILD/log(2), 'm-+', label='SB ILD')
#xlabel('P/N_0 (dB)')
#ylabel('Information Rate (bits/sec)')
#title('Channel ' + channel_index + ', Scheme ' + str(scheme_index))
#legend()
#
#show()

#%%
#
#DATA_DICT = {};
#DATA_DICT['SIM'] = SIM;
#DATA_DICT['SNRdB'] = SNRdB_vec;
#DATA_DICT['CHANNEL'] = CH;
#DATA_DICT['SCHEME'] = {};
#DATA_DICT['SCHEME']['PARAMS'] = SCHEMES_DICT[scheme_index];
#DATA_DICT['SCHEME']['RATE'] = {}
#DATA_DICT['SCHEME']['RATE']['R_OPTRX'] = RATE_OPT;
#DATA_DICT['SCHEME']['RATE']['R_SB_JLD'] = RATE_SB_JLD;
##DATA_DICT['SCHEME']['RATE']['R_SB_ILD1'] = R_SB_ILD1;
##DATA_DICT['SCHEME']['RATE']['R_SB_ILD2'] = R_SB_ILD2;

#%%

#def plot_spectrum(H_TX, H_CH, SCH, CH, SIM):
#[Sigma_X_NORMALIZED, layer] = power_alloc(H_TX, SCH, SIM);
#K_prime = length(layer);
#SELECTED_VECS = zeros(1,K_prime);
#for k in range(K_prime):
#    SELECTED_VECS(k) = layer{k}(2);
#
##%%
#N_FFT = SIM.T_SIMULATION/SIM.dt;
#f = linspace(-SIM.F_samp/2,SIM.F_samp/2,N_FFT);
#
#H_TX_SELECTED = H_TX(:,SELECTED_VECS);
#FT_IN = fftshift( fft(H_TX_SELECTED)*SIM.dt ,1);
#
#line_color = {'b',[0 0.5 0],'r','m'};
#
#W_base = SCH(1); a_base = SCH(2); K_prime = SCH(3); fc_base = SCH(4);
#f_max = (fc_base+W_base/2)*a_base^(K_prime-1);
#F = f_max * max(CH.alpha);
#
#H_CHTX = H_CH * H_TX_SELECTED;
#FT_OUT = fftshift( fft(H_CHTX)*SIM.dt , 1);
#
#figure(20)
#clf(20)
#subplot(K_prime+1,1,1)
#for k=1:K_prime
#    plot(f,abs(FT_IN(:,k)),'Color',line_color{k},'LineWidth',2)
#    hold on
#    axis([-F F 0 1])
#end
#hold off
#xlabel('f')
#ylabel('|FT|')
#title('Input Spectrum')
#
#for k=1:K_prime:
#    subplot(K_prime+1,1,k+1)
#    plot(f,abs(FT_OUT(:,k)),'Color',line_color{k},'LineWidth',2)
#    axis([-F F 0 1])
#    xlabel('f')
#    ylabel('|FT|')
#    title(sprintf('Output Spectrum (k''=%d)',k))
#
#set(20,'Position',[2*480 50 2*480 2*470]);figure(20)