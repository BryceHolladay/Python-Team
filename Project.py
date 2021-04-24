## ME396P Keith Lane Final Project
## Team Members: Ali Bhaiwala, Sarah Hildreth, Bryce Holladay
## Spring 2021

from scipy import signal
from scipy import fft
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
The purpose of this code is to analyze ground reaction force (GRF) data from a split-belt treadmill to 
identify the mean peak GRF values from the vertical GRFs while ignoring any peaks that count as 
'cross-over steps,' or steps by the subject on the incorrect treadmill.
"""

# Function that takes in a csv file and outputs an average max peak, as well as saves a graph
def find_max(file):
    csv = pd.read_csv(file, header=[3], nrows=36000) #Header [3] to only include 4th row of column titles
    
    frame = csv['Frame'].tolist()   
    fzR = csv['Fz'].tolist() #creating a list of all Fz values on right foot
    fzL = csv['Fz.1'].tolist() #create list of all Fz values on left leg
    
    #Removes units element from each column
    frame.pop(0)
    fzR.pop(0)
    fzL.pop(0)
    
    #changing the type of fzR and fzL to a float
    for i in range(len(fzR)):
        fzR[i] = float(fzR[i])
        
    for i in range(len(fzL)):
        fzL[i] = float(fzL[i])
    
    def project_filter(array, order, cutoff, filt):
    
        """
        This function performs a butterworth filter. Input parameters are the data to be 
        filtered, the order of the filter, the cutoff frequency(ies), and the type of filter

        Type of filter options include:
            'highpass'- highpass filter
            'lowpass'- lowpass filter
            'bandpass'- bandpass filter (requires range of frequencies)
            'bandstop'- bandstop filter(requires range of frequencies)

            ** default is lowpass filter

        """

        if filt == '':
           filt = 'lowpass' 

        b, a = signal.butter(order, cutoff, filt)

        filtarray= signal.filtfilt(b, a, array)

        return filtarray


    def avg_peak_isolation(array):
        # This method takes any peaks within 25% of the mean of the top 10 peaks
        pks_locs= signal.find_peaks(array) #gives the indicies of the peaks in array

        pks_all=[]
        for i in pks_locs:
            pks_all.append(array[i]) #creates array of values of each peak in array

        pks_all.sort()
        pks_all.reverse()
        pks_top10=pks_all[0:10]
        top10_mean= np.mean(pks_top10)


        good_pks=[]
        for peak in pks_all:
            if peak >= (0.75*top10_mean) and peak <= (1.25*top10_mean):
                good_pks.append(peak)

        return good_pks

    def weight_cutoff_isolation(array, mass, mult):
        # This method takes any peaks that are greater than (mult) times the subject's body weight
        # mass should be in kilograms

        weight= 9.81*mass
        pks_locs= signal.find_peaks(array, mult*weight) #returns array of indicies of the peaks in array

        good_pks=[]
        for i in pks_locs:
            good_pks.append(array[i]) #creates array of all peak values in array greater than (mult)*BW

        return good_pks

    return

find_max('S4_T21.csv')
