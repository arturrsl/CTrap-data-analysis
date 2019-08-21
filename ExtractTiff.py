#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:18:11 2019

@author: Artur K.

"""
import os
import numpy as np
import sys
import matplotlib.pyplot as plt

from lumicks import pylake
#import cv2 as cv
#matplotlib inline




folder =  r'/Users/Artur/Desktop/data/20190416a'

filenames = os.listdir(folder)          # all files in the chosen folder
Filenames = []                          

for filename in filenames:              # selection of .H5 files only
    if filename[-3:] == '.h5':
        Filenames.append(filename)

"""
 # Extract kymograph from all files in one folder and save them as TIFF
"""

for i in range(len(Filenames)):
    print (Filenames[i])
   
    #Load Marker File containing both kymograph and force extension data
    name = str(Filenames[i])
    file = pylake.File(name)
    #file = pylake.File("20190517-192654 ATP Kymograph 22.h5")
    print("Kymo I.D. = " + str(list(file.kymos)))
    print("FD Curve I.D. = " + str(list(file.fdcurves)))
    #print(file) #textual representation of the contents of a file


   

    kymo_names = list(file.kymos)                #reference to Kymo I.D 
    kymo = file.kymos[kymo_names[0]]
    plt.figure(figsize=(40,10))                  #changes plot size
    #kymo.plot_rgb()                             
    kymo2 = kymo.plot_green()
    #print("Start timestamp = " + str(kymo.start) + " ns")
    #print("End timestamp = " + str(kymo.stop) + " ns")
    time = round((kymo.stop - kymo.start) /1000000000,2)
    #round(number[, ndigits])
    print("Imaging time = " + str(time) + " s")
    
    index = kymo_names[0]                        #first index from the list of kymographs within one file    
    name = 'Kymo' + index + '.tiff'

    kymo.save_tiff(str(name))
    #kymo.save_tiff("Kymo22.tiff")







