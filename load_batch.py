# -*- coding: utf-8 -*-
"""
Created on Thu May 07 15:27:37 2015

@author: Hassan

Compare plots
Reading a structure from pickle (Python)

"""
#%%

#from multilayer5 import load_data
from numpy import *


#%% Get a batch of files (matching certain filename pattern)
import os
# get current working directory
cwd = os.getcwd()
# get all files in a given directory (path)
path = cwd + '\\data'
#path = cwd
#path = '.'

def get_sim_file_list(sim_filename_pattern):
    file_list = os.listdir(path)
    # get just the simulation output files
    import re # regular expression (regex)
    #sim_filename_pattern = '[0-9]{6}_py$' # string that ends with: six digits followed by '_py'
    sim_file_list = [f for f in file_list if re.search(sim_filename_pattern,f)!=None] # find files with the specified name pattern

    return sim_file_list

#sim_filename_pattern = '[0-9]{6}_py$' # string that ends with: six digits followed by '_py'
sim_filename_pattern = 'sarajevo_py$'
#sim_filename_pattern = '.mat$' # string that ends with: .mat
sim_file_list = get_sim_file_list(sim_filename_pattern)

#%% Index the filenames by 'channel' and 'scheme'
def index_filenames(sim_file_list):
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

    return filenames

if __name__ == '__main__':    
    filenames = index_filenames(sim_file_list)

#%%
#CH_LABELS = range(1,5+1)
#CH_LABELS = list(sort(filenames.keys()))
CH_LABELS = ['A', 'B', 'C', 'D', 'E']
SCH_LABELS = range(1,3+1)
from itertools import product # Caterzian/Cartesian product
for ch_index, sch_index in product(CH_LABELS, SCH_LABELS):
    #print('Channel %s, Scheme %d' %((str(ch_index), sch_index)))
    pass
	
    
#%% Load results from files

def get_results(filenames,file_type='py'):
    runtimes = []
    RESULTS = {}
    for ch_idx in CH_LABELS:
        RESULTS[('CHANNEL',ch_idx)] = {}
        for sch_idx in SCH_LABELS:
            print("Loading %s" % filenames[ch_idx][sch_idx])
            # load the data from files
            if file_type=='py': # load the appropriate function for python (pickle)
                import multilayer5
                load_data = multilayer5.load_data
            if file_type=='mat': # load the appropriate loading function for mat files
                import mat_interface
                load_data = mat_interface.load_data
            DATA = load_data(path+'\\'+filenames[ch_idx][sch_idx])
    		
            
            # extract the information rates
            RATES = {}
            for RX in DATA['RX'].keys():
                RATES[RX] = DATA['RX'][RX]['RATE']
            
            # store results            
            RESULTS[('CHANNEL',ch_idx)][('SCHEME',sch_idx)] = RATES
            RESULTS['SNRdB'] = DATA['SNRdB']; # rewriting is ok, since we run all simulations for the same SNR values
            
            runtimes.append(DATA.get('RUNTIME')) # If no 'RUNTIME', *None* will be rerturned
    
    #return RESULTS, runtimes
    return RESULTS

if __name__ == '__main__':
    RESULTS = get_results(filenames)

#%%
def load_batch(sim_filename_pattern,file_type='py'):
    sim_file_list = get_sim_file_list(sim_filename_pattern)
    filenames = index_filenames(sim_file_list)
    RESULTS = get_results(filenames,file_type)
    return RESULTS
    
#%%

if __name__ == '__main__':
    SNRdB_vec = RESULTS['SNRdB']
    
    LineWidth = 2;
    MarkerSize = 10;
    FontSize = 17;
    FontSizeLegend = 14;
    
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

    
    CHANNEL_LIST = ['E']
    SCHEME_LIST = [1, 2, 3]
    
    for ch_idx in CHANNEL_LIST:
        for sch_idx in SCHEME_LIST:
            figure(1000+sch_idx)
            plot_rate(ch_idx,sch_idx,'OPT'    , 'Optimal' , linestyle_list[1][1])
            plot_rate(ch_idx,sch_idx,'EB'     , 'EB'      , linestyle_list[1][2])
            plot_rate(ch_idx,sch_idx,'SB_JLD' , 'SB JLD'  , linestyle_list[1][3])
            plot_rate(ch_idx,sch_idx,'SB_ILD' , 'SB ILD'  , linestyle_list[1][4])
            title('Channel ' + str(ch_idx) + ', Scheme ' + str(sch_idx))
            xlabel('SNR (dB)')
            ylabel('Information Rates (nats/sec)')
            legend(loc='lower right')
    
    show()
