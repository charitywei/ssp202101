import numpy as np
import pandas as pd
import scipy.signal as sg
from scipy.stats import kstest
import matplotlib.pyplot as plt
from moonpy import *

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

        start = int(input('Kernel Size Start:'))
        end   = int(input('Kernel Size End:'))
        step  = int(input('Kernel Size Step:'))
        variable_input = input('Plots? y/n:')
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
                if variable_input == 'y':
                    plt.show()
                else:
                    pass