#!/usr/bin/env python
# coding: utf-8

# #                                      C-Trap kymograph analysis
# 
# July 2019
# Artur Kaczmarczyk

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import sys
from lumicks import pylake

import peakutils
from peakutils.plot import plot as pplot

import cv2


import PIL 
from PIL import Image

get_ipython().run_line_magic('matplotlib', 'inline')


# # Loading the .h5 file and inspecting its content
# 
# 1. Make sure that this notebook is in the same folder as the file that is about to be analyzed. 
# 2. Copy the entire name of the file below.
# 
# The .H5 file exported from Bluelake should contain at least one kymograph that has an assigned number (Kymo I.D). The script will further automatically call the kymograph by this number.

# In[2]:


file = pylake.File("20190123-195538 Kymograph 11.h5")
print("Kymo I.D. = " + str(list(file.kymos)))

# Loading the Kymo I.D. as a variable so there is no need to type it manually every time in other windows
kymo_number = (list(file.kymos))   

# Inspecting the content of the .h5 file
print(file)


# In[ ]:





# # Plotting and saving the kymograph together with the force measurement 
# 
# Kymograph is loaded as .png by default. The image's brightness is improved (vmax=N) to visualize the green signal. Next, force detection is downscaled and co-plotted with the kymograph.
# 
# 1. Change "green" to "red" or "blue" depending on which dye was used in the experiment !!
# 2. Adjust the PNG brightness by correcting the "vmax" parameter (10. line)

# In[3]:


index = 0 # keep ZERO when there is only one kymograph within the .h5 file, increments +1 would correspond with list indexes for remaining kymographs generated in the same .h5 file
kymo = file.kymos[(kymo_number[index])]


# Loading kymo and correcting brightness of the green signal

#image_rgb = kymo.plot_rgb()
green_img = kymo.green_image 
#max_px = np.max(green_img)
#print(max_px)           

ax = plt.figure(figsize=(30,20))
kymo.plot_green(vmax=6)  # decrease the value of VMAX if the obtained PNG image is too dimmed

axes = plt.gca()
xmin, xmax = axes.get_xlim()
plt.savefig('Kymo' + (kymo_number[index]) + '.png', dpi=300, bbox_inches='tight')



# Co-plotting the force

forcex = file['Force HF']['Force 2x']    # force in the x direction (pN)

sample_rate = forcex.sample_rate         # downsampling the force, nanostage position and time
downsampled_rate = 100 # Hz
forcex_downsamp = forcex.downsampled_by(int(sample_rate/downsampled_rate))

time = forcex.timestamps/1e9             # time traces (seconds)
time = time - time[0]
time_downsamp = forcex_downsamp.timestamps/1e9
time_downsamp = time_downsamp - time_downsamp[0]

#forcex.plot(label="Original")           # raw force data without downsampling
forcex_downsamp.plot(color='r',label="force downsampled to 100 Hz")
plt.ylabel('Force (pN)')
plt.xlim([0,max(time)])
plt.legend()

plt.savefig('KymoForce_' + (kymo_number[index]) + '.png', dpi=300, bbox_inches='tight')


# In[4]:


# Saving the PNG file and raw TIFF file 

kymo.save_tiff('Kymo_' + (kymo_number[index]) + '.tiff')
#kymo.save_tiff('Kymo_' + (kymo_number[index]) + '_8bit.tiff', dtype=np.int8) #necessary for ridge detection analysis
#plt.figure(figsize=(30,20))
#plt.xlim(xmin,xmax)
#file.force2x.plot()


# # Generating a cummulative intensity profile of the (sliced) kymograph
# The intensity of signal in green and blue channel is scanned along vertical lines of the kymograph. The intensity profiles are subsequently added and normalized and cummulative intensity profile is plotted.  
# 
# 1. Select manually which part of the kymograph needs to be analyzed by typing a start and end point (in seconds). 
# 

# In[4]:


green_img = kymo.green_image

line_time = time_downsamp[-1]/green_img.shape[1]
print('Line time = ' + str(round(line_time*1000)) + ' ms.')


# Type manually the timepoints (in seconds) that should be extracted from the full kymograph 

start = [270]
end = [550]

# Activate the two lines below if no slicing is required
#start = [0]
#end =  int(round((max(time_downsamp))))-1  # subtracting 1 sec in case the rounding approximates to the +1 larger integer 

cut_left = int(start/line_time)
cut_right = int(end/line_time)


green_img_sum = green_img[:,cut_left:cut_right].sum(axis=1)          
#green_img_sum = green_img[:,:].sum(axis=1)                # no slicing
norm_green_img_sum = green_img_sum/np.amax(green_img_sum)

blue_img = kymo.blue_image
blue_img_sum = blue_img[:,cut_left:cut_right].sum(axis=1)   
#blue_img_sum = blue_img[:,:].sum(axis=1)                  # no slicing
norm_blue_img_sum = blue_img_sum/np.amax(blue_img_sum)     


# Plotting together green and blue channel, normalized to the highest value of the profile at given channel 

fig, ax = plt.subplots(figsize=(30, 10))
ax.plot(norm_green_img_sum, 'g',norm_blue_img_sum, 'b')
ax.set_xlabel("pixel", fontsize=15)
ax.set_ylabel("normalized frequency", fontsize=15)

# Plotting separately 
#plt.figure(figsize=(30,10))
#plt.plot(green_img_sum/np.amax(green_img_sum), 'g'); 
#plt.figure(figsize=(30,10))
#plt.plot(blue_img_sum/np.amax(blue_img_sum), 'b');


# # Finding the centers of the beads and translating pixels to a genomic position
# 
# The script works with the peak detection package (download 'peakutils' library).
# 
# 1. Type the necessary information to translate pixel position into a genomic position (pixel size, bead size, DNA contour length).
# 2. Optionally, change the thresholds in the line corresponding with peak detection.
# 3. Select "flip" betwen 0 or 1 in case you need to reverse the polarity of the stands

# In[5]:


# Experimental parameters
ps = 100     # nm (pixel size) 
bs = 4420    # nm (bead size)
ds = 16490   # nm tethered DNA length, lambda DNA: 16490 nm
offset = 10  # number of plotted kymograph lines before/after the first/second bead-center peak, respectively
flip = 0     # flipping x axis, zero means No, 1 means Yes

# Generating an array with pixel indexes where a peak in intensity was detected in blue channel. (Set different thresholds if necessary!)
indexes = peakutils.indexes(norm_blue_img_sum, thres=0.3, min_dist=5)
print('Bead centers in lines ' + str(indexes) + ' correspond with -2.24 um and 18.70 um')


# Zooming to the region in between two blue peaks 
roi_left = indexes[0]-offset    
roi_right = indexes[1]+offset   

norm_green_img_sum_cut = norm_green_img_sum[roi_left:roi_right]
norm_blue_img_sum_cut = norm_blue_img_sum[roi_left:roi_right]

#length = indexes[1]-indexes[0] + 2*offset
#contour = ds + bs # lambda DNA contour length and one diameter of the bead



# Linear interpolation of line indexes into length in nanometers and then into kilobases

xp = [roi_left,roi_right]
fp = [-0.5*bs-offset*ps, ds+0.5*bs+offset*ps]
lines = np.arange(roi_left, roi_right,1)


position= np.interp(lines, xp, fp)/0.34/1000

if flip == 1:
    position = position[::-1]    #Activate this line if you want to swap the x axis

# Alternative wasy: position is calculated with the left blue peak as a reference and the genomic position is calculated based on the pixel size
# lines2 = np.arange(0,length,1)
# position2 = (0 - ps*offset - 0.5*bs + lines2*ps)/0.34/1000     
        

fig, ax = plt.subplots(figsize=(30, 10))
ax.plot(position, norm_green_img_sum_cut, 'g',position, norm_blue_img_sum_cut, 'b', linewidth=3.0)
ax.tick_params(direction='out', length=6, width=2, labelsize= 25,
               grid_color='r', grid_alpha=0.5)
ax.legend(['filament position ','centers of the beads'], fontsize=25)
ax.set_xlabel("genomic position (kB)", fontsize=25)
ax.set_ylabel("normalized frequency", fontsize=25)

plt.savefig('KymoIntensity' + (kymo_number[index]) + '.png', dpi=300, bbox_inches='tight')


# Detect peaks in the green channel and saving the result as CSV file

indexes2 = peakutils.indexes(norm_green_img_sum, thres=0.10, min_dist=3)
position_green = np.interp(indexes2, xp, fp)/0.34/1000    # interpolation 
#position_green = position_green[::-1]

position_green2 =  position_green[position_green>0]       # eliminating the green peaks coming from the bead1
position_green3 =  position_green2[position_green2<48.5]  # eliminating the green peaks coming from the bead2

if flip == 1:
    position_green3 = 48.500 - position_green3
    
    
print('Peaks in the green channel are on positions ' + str(position_green3) + ' kB' )
np.savetxt('KymoPeaks_' + (kymo_number[index]) + '.csv', position_green3, delimiter=" ")    # exports the peak position
np.savetxt('GreenIntensity_' + (kymo_number[index]) + '.csv',np.c_[position, norm_green_img_sum_cut], delimiter=",") # exports the entire green trace


# # Ridge detection
# 
# Locating the position (and movement) of fluoresently labeled molecule on the kymograph. Edge detection function is used (necessary 'cv2' library) that takes second derivative of each line. An edge / a ridge gives a peak. By finding the index of the peak, one can extract the genomic position.
# 

# In[6]:


import cv2
import matplotlib.pyplot as plt
from skimage.feature import hessian_matrix, hessian_matrix_eigvals


src_path = 'Kymo_' + (kymo_number[index]) + '_16bit.tif'    # opening the previously generated 8bit tiff file

# Ridge detection functions found online

def detect_ridges(gray, sigma):
    H_elems = hessian_matrix(gray, sigma=sigma, order='rc')
    maxima_ridges, minima_ridges = hessian_matrix_eigvals(H_elems)
    return maxima_ridges, minima_ridges

def plot_images(*images):
    images = list(images)
    n = len(images)
    fig, ax = plt.subplots(nrows=n, sharey=True)
    for i, img in enumerate(images):
        ax[i].imshow(img, cmap='gray')
        ax[i].axis('off')
    plt.subplots_adjust(left=0.03, bottom=0.03, right=2.5, top=0.97)
    plt.show()

img = cv2.imread(src_path, 0) # 0 imports a grayscale
if img is None:
    raise(ValueError(f"Image didn\'t load. Check that '{src_path}' exists."))

    
    
# dx is the first derivative of the image, ddx is the second derivative. Minima indicate the ridges

dx, ddx = detect_ridges(img[roi_left:roi_right,cut_left:cut_right], sigma=2)
plot_images(img, dx, ddx)       


# Get the minimum values of each column i.e. along axis 0
minInColumns = np.amin(ddx, axis=0)
#print('Min value of every column: ', minInColumns)

# NORMALIZATION (since I divide by the minimal negative number, second derivative will become positive and a valley will become a peak)
ddx_norm = ddx/minInColumns

#fig, ax = plt.subplots(figsize=(30, 10))  # ploting the raw image, 1st derivative image and 2nd derivative image
#ax.plot(ddx_norm[:,2],'b')                # co-plotting the peaks from a selected (or all) line (can be switched off)




[a,b] = ddx.shape  #b is the number of lines


# ACTIVATE IN ORDER TO SEARCH RIDGES IN KYMOGRAPHS WITH ONLY ONE SLIDING MOLECULE 
# [the loops detects MINIMUM VALUE IN EACH LINE TO DETECT THE EDGE]

#c = []
#for i in range(1,b):
    #print(np.argmin(ddx[:,i]))
#    c = np.append(c, np.argmin(ddx[:,i]))  #indexes of the lowest values in the column



# 
# SEARCH RIDGES IN KYMOGRAPHS WIYH MORE THAN ONE SLIDING MOLECULE

ridges = []      # Ridges are basically the indexes of elements detected by the peak detector of second derivative profiles

ridges_first = peakutils.indexes(ddx_norm[:,0], thres=0.2, min_dist=5)   # first line that will define how many ridges will be found in the entire trace
c = np.zeros((len(ridges_first),1))         # generate a single column with zeros
c = np.insert(c, 0, ridges_first, axis=1)   # put the array with the found positions of the peak as the first column of a newly made array "c"


def find_nearest(array, value):             # function that will search the value in an array that is the nearest to a defined number
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]



for i in range(1,b):          # looping through each line of the kymograph
 
    ridges = peakutils.indexes(ddx_norm[:,i], thres=0.2, min_dist=5)
    #print(ridges)
    add = []
    
    for j in range(0,len(ridges_first)):      # looping through each element of the vector with ridges position
        #print(ridges_first[j])
        nearest = find_nearest(ridges, ridges_first[j])
        
        if abs(nearest - ridges_first[j]) > 8:         # avoiding spikes (not sure if it's a good approach)
            nearest = ridges_first[j]
      
        #print(type(nearest))
        add = np.append(add,[nearest],axis=0)
        ridges_first[j] = add[j]
        
    c = np.insert(c, i, add, axis=1)
    

# Each row corresponds with a time-trace of one of the detected sliding molecules

c = c[:,:-1] # discarding the last column with zeros that was generated in an empty array 



#print(c[:,0:19])    
#print(c.shape)


fig, ax = plt.subplots(figsize=(18, 1))   # ploting the raw image, 1st derivative image and 2nd derivative image
ax.plot(c[0,:],'r')
plt.gca().invert_yaxis()
print((c.shape))


# ### Calculating the diffusion rate
# 

# In[ ]:





# In[7]:


g = c[0,:]
pos = np.around(g*(ps/1000),decimals =2)   # from pixel to micrometer knowing the pixel size

b = len(pos)
b = 200
print(g[0:10])
MSD = []
D = []
V = []
temp = []
shift = []
csum = []
for i in range(1,b):                    # looping through each column
 
    #MSD = (pos[:,i] - pos[:,0])**2     # or should I substract always the first data point in time zero?
    f = (pos[i] - pos[i-1])**2
    h = np.abs(pos[i] - pos[i-1])
    MSD = np.append(MSD,[f], axis =0)  
    shift = np.append(shift,[h], axis =0) 
    


for k in range(1,b):                    # looping through each column and summing all previous values along rows
    
    ff = (MSD[0:k].sum(axis=0))
    
    hh = (shift[0:k].sum(axis=0))
    D = np.append(D,[ff],axis=0)
    csum = np.cumsum(D)/k 
    V = np.append(V,[hh],axis=0)
    

time = line_time*np.arange(b-1)         # time of the tracked 1D diffusion
diff_rate = (csum[-1]/(time[-1]))       # diffusion rate counted as the max MSD in total time 
velocity = (V[-1]/time[-1])*1000/0.34   # total "route" divided by time (micrometers converted into nucleotides -> 1 bp = 0.34 nm)
dv_dt = (shift/line_time)*1000/0.34;
p,q = np.polyfit(time, csum, 1)         # fitting the linear function

print('Diffusion rate = ' + str(p/2) + ' um^2/s')


# An "interface" to matplotlib.axes.Axes.hist() method
n, bins, patches = plt.hist(x=dv_dt, bins=20, color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('dv/dt')
#plt.text(23, 45, r'$\mu=15, b=3$')
#maxfreq = n.max()
# Set a clean upper y-axis limit.
#plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

    
fig, ax = plt.subplots(figsize=(20, 5)) 

for z in range (0,len(D)):
    ax.plot(time, D, linewidth=3.0)      

#pp = np.random.rand(13)
#rr = np.cumsum((p**(0.5))*pp)


ax.tick_params(direction='out', length=6, width=2, labelsize= 25, grid_color='r', grid_alpha=0.5)
ax.set_xlabel("time (s)", fontsize=25)
ax.set_ylabel("MSD (um^2)", fontsize=25)

note = ('D =  ' + str(np.around(p/2,decimals =2))+' um^2/s')
note2 = ('V =  ' + str(np.around(velocity,decimals =2))+' nt/s')

ax.text(0.03, 0.9, note, fontsize=20)
ax.text(0.03, 1.2, note2, fontsize=20)


plt.savefig('MSD_time' + (kymo_number[index]) + '.png', dpi=300, bbox_inches='tight')


# Compiling time and diffusion rate into a single array that will be exported
diff_final = D.reshape(-1,1)
print(diff_final.shape)
print(time.shape)
time = time.reshape(-1,1)
S = np.column_stack((time,diff_final[:,0]))

np.savetxt('Diffusion_' + (kymo_number[index]) + '.csv', S, delimiter=",")    # exports the peak position
    



# In[182]:





# In[183]:





# In[184]:





# In[ ]:




