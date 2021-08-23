import numpy as np
import pandas as pd
import scipy.signal as sg
from scipy.stats import kstest
import matplotlib.pyplot as plt
from moonpy import *

def butter_lowpass_filter(data, cutoff = 1/3):
    # sample rate, /day
    # desired cutoff frequency of the filter, /day
    nyq = int(input('What is the nyq value?'))
    # Nyquist Frequency:if the data points are 2 minutes = 0.00138889 days, 
    # then the sampling frequency is 720 Hz. nyq is then half that value. 
    normal_cutoff = cutoff / nyq
    b, a = sg.butter(2, normal_cutoff, btype='low', analog=False)
    y = sg.filtfilt(b, a, data)
    return y

targets = [input('Which Star?')]
telescopes = ['TESS', 'Kepler', 'K2']


for target in targets:
    for telescope in telescopes:
        objects = MoonpyLC(targetID=target,telescope=telescope,load_lc='y',clobber='n')
        detrended = objects.detrend(dmeth = 'cofiam')
        fluxes  = objects.fluxes
        times   = objects.times

        fig = plt.figure(figsize=(8,6), dpi=300)
        ax = fig.add_subplot(1,1,1)
        for i in range(0,len(fluxes)):
            print(i,len(times[i]),len(fluxes[i]))
            ax.scatter(times[i],fluxes[i],s = 0.1, color = 'blue')
        plt.title(target)
        ax.set_xlabel(r'time (BTJD)')
        ax.set_ylabel(r'flux (arb.)')
        ax.set_xlim(np.min(times[0]),np.max(times[0]))

        ks = kstest(ls_medfluxes, np.random.normal(loc=1.0, scale=np.std(ls_medfluxes), size=len(ls_medfluxes)), N=20, alternative='two-sided', mode='auto')
        print(ks)

        start = int(input('Kernel Size Start:'))
        end   = int(input('Kernel Size End:'))
        step  = int(input('Kernel Size Step:'))
        variable_input = input('Plots? y/n:')
        if variable_input == 'y':
            plt.show()
        else:
            pass

        for i in range(0,len(fluxes)):
            for kern in np.arange(start, end, step):
                fig = plt.figure(figsize=(8,6), dpi=300)
                ax = fig.add_subplot(1,1,1)
                print(i,len(times[i]),len(fluxes[i]))
                med_fluxes = sg.medfilt(fluxes[i], kernel_size = kern)
                ax.set_xlim(np.min(times[i]),np.max(times[i]))
                ax.plot(times[i], med_fluxes, linestyle = '-', color = 'red')
                ax.set_title(kern)
                #fig.savefig('/Users/charitywei/Desktop/test01/' + target + '_'+ telescope + 'kernel size' + str(kern) + '_' + str(i) + 'medfilt.png', dpi=300)
        
        for i in range(0,len(fluxses)):
            fig = plt.figure(figsize=(8,6), dpi=300)
            ax = fig.add_subplot(1,1,1)
            ax.plot(kern,ks)
            plt.show()