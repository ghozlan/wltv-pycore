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
channel_index = 'B'

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

#%%

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
    
    figure()
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
    
    show()

#set(20,'Position',[2*480 50 2*480 2*470]);figure(20)

plot_spectrum(H_TX, H_CH, SCH, CH, SIM_DICT)