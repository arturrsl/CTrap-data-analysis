#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from lumicks import pylake


# ## Load FD curve

# In[ ]:


file = pylake.File("20190709-153836 pFM008 in NTM FD Curve 78.h5")
print("Kymo I.D. = " + str(list(file.fdcurves)))

# Loading the I.D. as a variable so there is no need to type it manually every time in other windows
curve_number = (list(file.fdcurves))   

# Inspecting the content of the .h5 file
#print(file)


# ## WLC - 7 parameter model
# This function returns the forces computed from a 7 parameter model of the WLC using the model by Bouchiat et al. Biophys J 76:409 (1999)
# 
# 

# In[ ]:


def wlc_7param(ext, Lp, Lc):
 
    # Lp = persistence length  (in micro-m)
    # Lc = contour length      (in micro-m)
    # T  = absolute temperature (in Kelvin)
    # k_B T in units of pN micro-m
    
    #Lc = 16.49
    T = 273
    Fwlc = []
    kT = 1.3806503 * 10**(-23) *T/(10**(-18)) 
    z_L = ext/Lc

    
    my_list = [1, -0.5164228, -2.737418, 16.07497, -38.87607,  39.49944, -14.17718]   # parameters from the paper

    a = np.asarray(my_list)
    Fwlc = 1/(4*(1 - z_L)**2)  - 1/4
    
    for i in range(0,len(a)):
        Fwlc = Fwlc + a[i] * z_L**[i]


    return Fwlc * kT/Lp;


# ## Plot the data and fit the persistence length

# In[ ]:


index = 0 # keep ZERO when there is only one kymograph within the .h5 file, increments +1 would correspond with list indexes for remaining kymographs generated in the same .h5 file
fd = file.fdcurves[(curve_number[index])]

force = fd.f
distance = fd.d


# Plot the stretching curve
plt.figure(figsize=(10, 7), dpi=70)
plt.scatter(distance.data[0:np.argmax(distance.data)], force.data[0:np.argmax(distance.data)], s = 12)
# Plot the refolding curve
#plt.scatter(distance.data[np.argmax(force.data):len(force.data)], force.data[np.argmax(force.data):len(force.data)], s = 12, color = 'silver')


# Choose the part of the data that will be fitted   
force_treshold = 10    # pN   
i = 0
while (force.data[i] < force_treshold):
    arg = force.data[index]
    i = i+1
print(i)     
input_data = distance.data[0:i]

# Least square method to fit the optimal persistence length and contour length 

popt, pcov = curve_fit(wlc_7param, input_data, force.data[0:i], bounds=([0.005,6.301],[0.060,6.801]))
perr = np.sqrt(np.diag(pcov)) #standard deviation errors on the parameters
print(perr[0]*1000)
print(perr[1])
plt.plot(input_data, wlc_7param(input_data, *popt), 'r-',linewidth=2.5, label='fit: a=%5.3f, b=%5.3f' % tuple(popt))


# Plot WLC with a given persistence length

#length = np.arange(3000, 16200, 45)/1000       #render an array of extension for dsDNA
#force_Calc = wlc_7param(length,0.045,16.49)     
#plt.plot(length, force_Calc,'g:', linewidth=3)
   

plt.xlabel('extension (\u03BCm)')    
plt.ylabel('force (pN)')   
plt.xlim(3, 7)
plt.ylim(-0.2, 30)

plt.rc('font', size=15)

leg = plt.legend(('fit WLC, PL = '+ str(round(popt[0]*1000))+ ' nm, CL = '+str(round(popt[1],2))+' \u03BCm','stretching curve'))
leg.get_frame().set_linewidth(0.0)

plt.savefig('fit_FD_' + (curve_number[index]) + '.png', dpi=300, bbox_inches='tight')


# In[ ]:




