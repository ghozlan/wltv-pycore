# -*- coding: utf-8 -*-
"""
Created on Thu May 07 15:27:37 2015

@author: Hassan

Compare plots
Reading a structure from pickle (Python)

"""
#%%

from multilayer4 import load_data
from numpy import *

#DD = load_data('channela_scheme1_125123_py')
#D = load_data('multilayer\\channele_scheme4_160436_py') # new

#%%

import os
# get current working directory
cwd = os.getcwd()
# get all files in a given directory (path)
path = cwd + '\\multilayer'
path = cwd
file_list = os.listdir(path)
# get just the simulation output files
#sim_filename_pattern = '[0-9]{6}_py$' # string that end withs: six digits followed by '_py'
sim_file_list = [f for f in file_list if re.search(sim_filename_pattern,f)!=None] # find files with the specified name pattern
#sim_file_list = [f for f in file_list if f.find('_py')!=-1] # does not work properly


#%%
filenames = {}
for filename in sim_file_list:
    # extract the channel and the scheme from the filename    
    filename_ch, filename_sch = filename.split('_')[0:2];
    filename_ch = filename_ch[len('channel'):].upper()
    filename_sch = int( filename_sch[len('scheme'):] )
    print('File  %s : Channel %s, Scheme %s' %( (filename,str(filename_ch),filename_sch) ))
    
    # store filenames in a dict indexed by channel and scheme    
    if filenames.get(filename_ch)==None: filenames[filename_ch] = {}
    filenames[filename_ch][filename_sch] = filename

#%%
#CH_LABELS = range(1,5+1)
#CH_LABELS = list(sort(filenames.keys()))
CH_LABELS = ['A', 'B', 'C', 'D', 'E']
SCH_LABELS = range(1,3+1)
from itertools import product # Caterzian/Cartesian product
for ch_index, sch_index in product(CH_LABELS, SCH_LABELS):
    print('Channel %s, Scheme %d' %((str(ch_index), sch_index)))
    
#%%
runtimes = []
RESULTS = {}
for ch_idx in CH_LABELS:
    RESULTS[('CHANNEL',ch_idx)] = {}
    for sch_idx in SCH_LABELS:
        print("loading %s" % filenames[ch_idx][sch_idx])
        # load the data from files
        DATA = load_data(filenames[ch_idx][sch_idx])
        
        # extract the information rates
        RATES = {}
        for RX in DATA['RX'].keys():
            RATES[RX] = DATA['RX'][RX]['RATE']
            
        RESULTS[('CHANNEL',ch_idx)][('SCHEME',sch_idx)] = RATES
        SNRdB_vec = DATA['SNRdB'];
        
        runtimes.append(DATA['RUNTIME'])

RESULTS_PY = RESULTS

#%%

LineWidth = 2;
MarkerSize = 10;
FontSize = 17;
FontSizeLegend = 14;

CHANNEL_LABEL = [' ', 'A', 'B', 'C', 'D', 'E']
from matplotlib.pyplot import *

def plot_rate(ch_idx,sch_idx,rate, label_str, linestyle_str='-'):
    line, = plot(SNRdB_vec,RESULTS[('CHANNEL',ch_idx)][('SCHEME',sch_idx)][rate], \
        linestyle_str, \
        linewidth=LineWidth, \
        markersize=MarkerSize, \
        label=label_str)
    #if style!=None: line.set_linestyle(linestyle_str);
        
#    plot(SNRdB_vec,RESULTS[('CHANNEL',ch_idx)][('SCHEME',sch_idx)]['R_OPT'], \
#        line_style[sch_idx][1],'LineWidth',LineWidth,'MarkerSize',MarkerSize)

linestyle_list = {}
linestyle_list[1] = [ '-', 'g*-', 'bs-', 'm+-', 'rx-', 'bo-' ];
linestyle_list[2] = [ '-', 'g*--', 'bs--', 'm+--', 'rx--', 'bo--' ];
linestyle_list[3] = [ '-', 'g*:', 'bs:', 'm+:', 'rx-', 'bo:' ];

ch_idx = 'E';
for sch_idx in range(1,3+1):
    figure(1000+sch_idx)
    plot_rate(ch_idx,sch_idx,'OPT'    , 'Optimal' , linestyle_list[1][1])
    plot_rate(ch_idx,sch_idx,'EB'     , 'EB'      , linestyle_list[1][2])
    plot_rate(ch_idx,sch_idx,'SB_JLD' , 'SB JLD'  , linestyle_list[1][3])
    plot_rate(ch_idx,sch_idx,'SB_ILD' , 'SB ILD'  , linestyle_list[1][4])
    legend(loc='lower right')
    title('Channel ' + str(ch_idx) + ', Scheme ' + str(sch_idx))
    xlabel('SNR (dB)')
    ylabel('Information Rates (nats/sec)')

show()

#import sys; sys.exit()
#%%
figure(77)
#plot(SNRdB_vec, RESULTS[('CHANNEL',ch_idx)][('SCHEME',sch_idx)]['R_OPT'])
line, = plot(SNRdB_vec, RESULTS[('CHANNEL',ch_idx)][('SCHEME',sch_idx)]['OPT'],linewidth=LineWidth,markersize=MarkerSize)
#title('Channel ' + str(CHANNEL_LABEL[ch_idx]) + ', Scheme ' + str(sch_idx))
xlabel('SNR (dB)')
ylabel('Information Rates (nats/sec)')
line.set_label('Opt RX')
legend()
##handles = get(1000+sch_idx,'Children'); set(handles(1),'FontSize',FontSizeLegend);    
#fig_filename = 'RATE_CH' + str(CHANNEL_LABEL[ch_idx]) + '_SCHEME' + str(sch_idx) + '.eps'
##print('-depsc', fig_filename)

