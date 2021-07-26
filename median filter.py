#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 15:39:11 2021

@author: charitywei
"""

import numpy as np
import pandas as pd
import scipy.signal as sg
from scipy.stats import kstest
import matplotlib.pyplot as plt
from moonpy import *

targets = ['HD 39060']
telescopes = ['TESS', 'Kepler', 'K2']


def butter_lowpass_filter(data, kernel_size = 999,cutoff = 1/5):
    # sample rate, /day
    # desired cutoff frequency of the filter, /day
    nyq = 0.5 * kernel_size  # Nyquist Frequency
    normal_cutoff = cutoff / nyq
    b, a = sg.butter(2, normal_cutoff, btype='low', analog=False)
    y = sg.filtfilt(b, a, data)
    return y

#detrending function
#def detrender(values):
    
    #detrended = []
    
    #Loop over the observation blocks
    #for i in range(0,len(values)):
        #Initial pass over the data block, exclude 3-sigma excursions from 
        #median calculation
        #valstd = np.std(values[i])
        #valmed = np.median(values[i])
        #a = np.where((abs(values[i] - valmed) / valstd) < 3.)
        
        #Calculate median on filtered data
        #valmed = np.median(values[i][a])
        
        #valfilt = np.asarray(values[i]/valmed)
        
        #detrended.append(valfilt)
    
    #detrended = np.asarray(detrended,dtype=object)
    
    #return detrended


for target in targets:
    for telescope in telescopes:
        try:
            objects = MoonpyLC(targetID=target,telescope=telescope,load_lc='y',clobber='n')
            detrended = objects.detrend(dmeth = 'cofiam')
            fluxes  = objects.fluxes
            times   = objects.times
            #detrend_fluxes = detrender(fluxes)
            
            #Make 7 rows into 1 row.
            #"+= list" is to connect items. 
            ls = list()
            for i in range(0,len(fluxes)):
                ls += list(fluxes[i])
            ls = np.array(ls)
            print(ls.shape)
            
            
            
            #fig = plt.figure(figsize=(8,6), dpi=300)
            #ax = fig.add_subplot(1, 1, 1)
            #ax.hist(ls, bins = 1000)
            #ax.hist(np.random.normal(loc=1.0, scale=np.std(ls), size=len(ls)), bins = np.linspace(0.95,1.05,100))
            #plt.show()
            fig = plt.figure(figsize=(8,6), dpi=300)
            ax = fig.add_subplot(1,1,1)
            for i in range(0,len(fluxes)):
                print(i,len(times[i]),len(fluxes[i]))
                ax.scatter(times[i],fluxes[i],s = 0.1, color = 'blue')
            plt.title(target)
            ax.set_xlabel(r'time (BTJD)')
            ax.set_ylabel(r'flux (arb.)')
            ax.set_xlim(1410, 1437)
            ax.set_ylim(5.3e6,5.45e6)
            #ax.set_xlim(np.min(times[0]),np.max(times[-1]))
            
            #Make 7 rows into 1 row.
            #"+= list" is to connect items. 
            #ls_medfluxes = list()
            ls_lpffluxes = list()
            
            for i in range(0,len(fluxes)):
                print(i,len(times[i]),len(fluxes[i]))
                #med_fluxes = sg.medfilt(fluxes[i], kernel_size = [999])
                lpf_fluxes = butter_lowpass_filter(fluxes[i])
                #ls_medfluxes += list(med_fluxes)
                ls_lpffluxes += list(lpf_fluxes)
                #ax.plot(times[i], med_fluxes, linestyle = '-', color = 'red')
                #ax.scatter(times[i], lpf_fluxes, s = 0.1, color = 'red')
                #ax.plot(times[i], lpf_fluxes, linestyle = '-', color = 'red')
            #ls_medfluxes = np.array(ls_medfluxes)
            ls_lpffluxes = np.array(ls_lpffluxes)
            #print(ls_lpffluxes.shape)
            #print(ls_medfluxes.shape)
            #ax.set_xlim(1410, 1438)
            #ax.set_ylim(-1e-10, 1e-10)
            #plt.savefig('/Users/charitywei/Desktop/HD39060_lpf.png')
            #fig.savefig('/Users/charitywei/Desktop/' + target + '_'+ telescope + '_'+ 'medfilt.png', dpi=300)
            #plt.show()
            
            ls_times = list()
            for i in range(0,len(times)):
                ls_times += list(times[i])
            ls_times = np.array(ls_times)
            print(ls_times.shape)
            
            fig = plt.figure(figsize=(8,6), dpi=300)
            ax = fig.add_subplot(1,1,1)
            #residuals = ls/ls_medfluxes
            residuals = ls/ls_lpffluxes
            ax.set_xlim(1410,1437)
            ax.set_ylim(0.9,1.1)
            #ax.set_title(target + '_' + 'Median Filtered')
            ax.set_title(target + '_' + 'Low Pass Frequency Filtered')
            ax.plot(ls_times, residuals, linestyle = '-', color = 'red')
            plt.show()
            
            fig, [ax1, ax2] = plt.subplots(2, sharex=True, dpi=300)
            ax1.scatter(ls_times, ls, s = 0.1)
            ax2.plot(ls_times, residuals, linestyle = '-')
            ax1.set_xlim([1410, 1437])
            ax1.set_ylim([5.32e6, 5.42e6])
            ax1.set_title(target + '_' + 'Raw Data')
            ax2.set_xlim([1410, 1437])
            ax2.set_title(target + '_' + 'Low Pass Frequency Filtered Residuals')
            #ax2.set_title(target + '_' + 'Median Filtered Residuals')
            plt.show()
            
            fig = plt.figure(figsize=(8,6), dpi=300)
            ax = fig.add_subplot(1, 1, 1)
            ax.hist(residuals, bins = np.linspace(0.99,1.01,1000), alpha = 0.5)
            ax.hist(np.random.normal(loc=1.0, scale=np.std(residuals), size=len(residuals)), bins = np.linspace(0.99,1.01,1000), alpha = 0.5)
            ax.set_title(target + '_' + 'Low Pass Frequency Filtered Residuals Histogram')
            #ax.set_title(target + '_' + 'Median Filtered Residuals Histogram')
            ax.set_xlim(0.995,1.005)
            plt.show()
            
            
        
            ks = kstest(ls_medfluxes, np.random.normal(loc=1.0, scale=np.std(ls_medfluxes), size=len(ls_medfluxes)), N=20, alternative='two-sided', mode='auto')
            print(ks)
        except:
            print('No data available for '+target+' in '+telescope+'.')

