#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 15:37:03 2021

@author: charitywei
"""

import cv2
import numpy as np
import glob
import os

dir_name = '/Users/charitywei/Desktop/asiaa/animation/KIC/lpf subplots/'
# Get list of all files in a given directory sorted by name
list_of_files = sorted( filter( os.path.isfile,
                        glob.glob(dir_name + '*') ) )
# Iterate over sorted list of files and print the file paths 
# one by one.

imgs_path = []
for file_path in list_of_files:
    #print(file_path) 
    imgs_path.append(file_path)
    
img_array = []
#for filename in glob.glob(imgs_path):
for filename in imgs_path:
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)
 
 
out = cv2.VideoWriter('KIC_lpf_subplots.mp4',cv2.VideoWriter_fourcc(*'DIVX'), 1, size)
 
for i in range(len(img_array)):
    out.write(img_array[i])
out.release()