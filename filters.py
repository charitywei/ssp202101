#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 20:01:22 2021

@author: charitywei
"""

import numpy as np
import pandas as pd
import scipy.signal as sg
from scipy import ndimage, misc
from scipy.stats import kstest
from scipy.interpolate import interp1d
from scipy import stats
from fractions import Fraction
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from moonpy import *

targets = [input('Which target?')]
telescopes = ['TESS', 'Kepler', 'K2']

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro')

def init():
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(-1, 1)
    return ln,

def update(frame):
    xdata.append(frame)
    ydata.append(np.sin(frame))
    ln.set_data(xdata, ydata)
    return ln,

def butter_lowpass_filter(data, cutoff):
    # sample rate, /day
    # desired cutoff frequency of the filter, /day
    nyq = 360  # Nyquist Frequency
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

savepath = '/Users/charitywei/Desktop/asiaa/animation/'

for target in targets:
    for telescope in telescopes:
        try:
            objects = MoonpyLC(targetID=target,telescope=telescope,load_lc='y',clobber='n')
            #detrended = objects.detrend(dmeth = 'cofiam')
            fluxes  = objects.fluxes
            
            times   = objects.times
            errors  = objects.errors
            norm_errors = errors/fluxes
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
            ax.set_xlim([np.min(times[0]),np.max(times[0])])
            ax.set_ylim([np.nanmin(fluxes[0]), np.nanmax(fluxes[0])])
            print('work line 92')
            
            #remove outliers
            for i in range(0, len(fluxes)):
                median = ndimage.median_filter(fluxes[i], size = 19)
                fluxes_std = np.nanstd(fluxes[i])
                diff = (fluxes[i] - median)/fluxes_std
                print('work line 96')
                outlier_idxs = np.where(abs(diff) >= 3)[0]
                inliers_idxs = np.where(abs(diff) < 3)[0]
                plt.scatter(times[i],fluxes[i],color='blue')
                plt.scatter(times[i][outlier_idxs],fluxes[i][outlier_idxs],color='red')
                plt.show()
                fluxes[i][outlier_idxs] = np.nan
                print('work line 101')
                
                #outliers = np.where(abs(fluxes[i] - np.median(fluxes[i])) / np.std(fluxes[i]))
                
                
                
            #sigma clip
            #for i in range(0, len(fluxes)):
                #fluxes[i] = stats.sigmaclip(fluxes[i], low=10.0, high=10.0)
            
            #plot the masked values (np.nans) on top of the original time series in a different colour. 
            
            
            for i in range(0,len(fluxes)):
                ls_pvalue = list()
                ls_cutoff = list()
                
                cutoff = np.arange(5, 200, 25)
                kernel_size = np.arange(55,555,50)
                for kern in kernel_size:
                #for lpf_cut in cutoff:
                    
                    print(i,len(times[i]),len(fluxes[i]))
                    #lpf_fluxes = butter_lowpass_filter(fluxes[i], cutoff = lpf_cut)
                    #lpf_residuals = fluxes[i]/lpf_fluxes
                    
                    print('work line 115')
                    notnan = np.isfinite(fluxes[i])
                    med_fluxes = ndimage.median_filter(fluxes[i][notnan], size = kern)
                    func = interp1d(times[i][notnan],med_fluxes)
                    actual_med_fluxes = func(times[i])
                    med_residuals = fluxes[i]/actual_med_fluxes
                    print('work line 118')
                    #med_fluxes = sg.medfilt(fluxes[i], kernel_size = kern)
                    #med_residuals = fluxes[i]/med_fluxes
                    #ks = kstest(lpf_residuals, np.random.normal(loc=1.0, scale=np.median(norm_errors[i]), size=len(lpf_residuals)), N=20, alternative='two-sided', mode='auto')
                    ks = kstest(med_residuals, np.random.normal(loc=1.0, scale=np.median(norm_errors[i]), size=len(med_residuals)), N=20, alternative='two-sided', mode='auto')
                    print(ks.pvalue)
                    print('work line 124')
                    #ls_pvalue.append(ks.pvalue)
                    #ls_cutoff.append(lpf_cut)
                    #ax.plot(times[i],lpf_residuals)
                    fig, [ax1, ax2] = plt.subplots(2, sharex=True, dpi=300)
                    ax1.scatter(times[i],fluxes[i], s = 0.1)
                    ax1.plot(times[i], actual_med_fluxes, linestyle = '-', color = 'red')
                    ax1.set_ylabel('Fluxes(arb.)')
                    ax1.set_title('Median Filtered' + '(' + 'Sec' + str(i) + ')(' + 'Kernel Size:' + str(kern) + ')')
                    ax2.scatter(times[i], med_residuals, s = 0.1)
                    ax2.set_xlabel('Time(BTJD)')
                    #ax2.plot(times[i], med_med_res, linestyle = '-', color = 'red')
                    fig.savefig(savepath + target + 'Sec' + str(i) + 'Kern' + str(kern) + 'Medfilt_subplots.png', dpi=300)
                    plt.show()
                    print('work line 130')
                    fig = plt.figure(figsize=(8,6), dpi=300)
                    ax = fig.add_subplot(1, 1, 1)
                    ax.hist(med_residuals, bins = np.linspace(0.99,1.01,1000), alpha = 0.5)
                    ax.hist(np.random.normal(loc=1.0, scale=np.std(fluxes[i][notnan]/np.median(fluxes[i][notnan])), size=len(fluxes[i][notnan])), bins = np.linspace(0.99,1.01,1000), alpha = 0.5)
                    ax.set_title('Median Filtered Residuals Histogram' + '(' + 'Sec' + str(i) + ')(' + 'Kernel Size:' + str(kern) + ')')
                    ax.set_xlim(1-3*np.std(fluxes[i][notnan]/np.median(fluxes[i][notnan])), 1+3*np.std(fluxes[i][notnan]/np.median(fluxes[i][notnan])))
                    fig.savefig(savepath + target + 'Sec' + str(i) + 'Kern' + str(kern) +'medfilt_hist.png', dpi=300)
                    plt.show()
                    #keep_going = input('Do you want to keep going?')
                    #if keep_going == 'n':
                        #raise Exception()
                #ax.set_xlim(np.min(times[i]),np.max(times[i]))
                #ax.plot(ls_cutoff, ls_pvalue, linestyle = '-', color = 'red')
                #ax.plot(kern, ls_pvalue, linestyle = '-', color = 'red')
                #plt.show()

                
                for lpf_cut in cutoff:
                    notnan = np.isfinite(fluxes[i])
                    lpf_fluxes = butter_lowpass_filter(fluxes[i][notnan], cutoff = lpf_cut)
                    func = interp1d(times[i][notnan],lpf_fluxes)
                    actual_lpf_fluxes = func(times[i])
                    lpf_residuals = fluxes[i]/actual_lpf_fluxes
                    fig, [ax1, ax2] = plt.subplots(2, sharex=True, dpi=300)
                    ax1.scatter(times[i],fluxes[i], s = 0.1)
                    ax1.plot(times[i], actual_lpf_fluxes, linestyle = '-', color = 'red')
                    ax1.set_ylabel('Fluxes(arb.)')
                    ax1.set_title('Low Pass Frequency Filtered' + '(' + 'Sec' + str(i) + ')(' + 'Cutoff:' + str(lpf_cut) + ')')
                    ax2.scatter(times[i], lpf_residuals, s = 0.1)
                    ax2.set_xlabel('Time(BTJD)')
                    fig.savefig(savepath + target + 'Sec' + str(i) + 'Cutoff' + str(lpf_cut) +'lpf_subplots.png', dpi=300)
                    plt.show()
                    fig = plt.figure(figsize=(8,6), dpi=300)
                    ax = fig.add_subplot(1, 1, 1)
                    ax.hist(lpf_residuals, bins = np.linspace(0.99,1.01,1000), alpha = 0.5)
                    ax.hist(np.random.normal(loc=1.0, scale=np.std(fluxes[i][notnan]/np.median(fluxes[i][notnan])), size=len(fluxes[i][notnan])), bins = np.linspace(0.99,1.01,1000), alpha = 0.5)
                    ax.set_title('Low Pass Frequency Filtered Residuals Histogram' + '(' + 'Sec' + str(i) + ')(' + 'Cutoff:' + str(lpf_cut) + ')')
                    ax.set_xlim(1-3*np.std(fluxes[i][notnan]/np.median(fluxes[i][notnan])), 1+3*np.std(fluxes[i][notnan]/np.median(fluxes[i][notnan])))
                    fig.savefig(savepath + target + 'Sec' + str(i) + 'Cutoff' + str(lpf_cut) +'lpf_hist.png', dpi=300)
                    plt.show()
        except:
            print('No data available for '+target+' in '+telescope+'.')