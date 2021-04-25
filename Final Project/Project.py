## ME396P Keith Lane Final Project
## Team Members: Ali Bhaiwala, Sarah Hildreth, Bryce Holladay
## Spring 2021

from scipy import signal
from scipy import fft
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import os

"""
The purpose of this code is to analyze ground reaction force (GRF) data from a split-belt treadmill to 
identify the mean peak GRF values from the vertical GRFs while ignoring any peaks that count as 
'cross-over steps,' or steps by the subject on the incorrect treadmill.
"""


def butter_filter(array, order, cutoff, dt, filt='lowpass'):
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

    cutoff_new = cutoff / (math.sqrt(2) - 1) ** (0.5 / order)

    b, a = signal.butter(order, 2 * dt * cutoff, filt)

    filtarray = signal.filtfilt(b, a, array)

    return filtarray


def avg_peak_isolation(array):
    # This method takes any peaks within 25% of the mean of the top 10 peaks
    pks = signal.find_peaks(array)  # gives the indicies of the peaks in array
    pks_locs = pks[0]
    pks_all = []

    for i in pks_locs:
        pks_all.append(array[i])  # creates array of values of each peak in array

    pks_all.sort()
    pks_all.reverse()
    pks_top10 = pks_all[0:10]
    top10_mean = np.mean(pks_top10)

    good_pks_top10 = []
    for peak in pks_all:
        if peak >= (0.75 * top10_mean) and peak <= (1.25 * top10_mean):
            good_pks_top10.append(round(peak, 3))

    return good_pks_top10


def weight_cutoff_isolation(array, mass, mult=2):
    # This method takes any peaks that are greater than (mult) times the subject's body weight
    # Returns a list of floats
    # mass is in kilograms

    weight = 9.81 * mass
    pks = signal.find_peaks(array, mult * weight)  # returns array of indicies of the peaks in array
    pks_locs = pks[0]
    good_pks_BWcut = []
    for i in pks_locs:
        good_pks_BWcut.append(round(array[i], 3))  # creates array of all peak values in array greater than (mult)*BW

    return good_pks_BWcut


# Function that takes in a csv file and subject's mass (Kg) and outputs an average max peak, as well as saves a graph
# noinspection PyUnresolvedReferences
def find_max(file, subjectmass):
    csv = pd.read_csv(file,
                      names=['fzL', 'fzR'],
                      usecols=[4, 13],
                      header=0,
                      skiprows=5,
                      dtype={'fzL': np.float64, 'fzR': np.float64},
                      nrows=35999)

    fzR = csv['fzR'].tolist()  # creating a list of all Fz values on right foot
    fzL = csv['fzL'].tolist()  # create list of all Fz values on left leg

    fzR = np.negative(fzR)
    fzL = np.negative(fzL)

    # filter the left and right GRFs with butterworth filter
    fzR_filt = butter_filter(fzR, 4, 15, 1 / 960, 'lowpass')
    fzL_filt = butter_filter(fzL, 4, 15, 1 / 960, 'lowpass')

    fzR_top10 = avg_peak_isolation(fzR_filt)  # all peaks accepted by avg peak cutoff method
    fzR_BWcut = weight_cutoff_isolation(fzR_filt, subjectmass)  # all peaks accepted by BW cutoff method
    fzL_top10 = avg_peak_isolation(fzL_filt)  # same but for left foot
    fzL_BWcut = weight_cutoff_isolation(fzL_filt, subjectmass)  # same but for left foot

    fzR_top10_mean = round(np.mean(fzR_top10), 3)  # mean peak value for top 10 avg method (right)
    fzR_BWcut_mean = round(np.mean(fzR_BWcut), 3)  # mean peak value for BW cutoff method (right)
    fzL_top10_mean = round(np.mean(fzL_top10), 3)  # mean peak value for top 10 avg method (left)
    fzL_BWcut_mean = round(np.mean(fzL_BWcut), 3)  # mean peak value for BW cutoff method (left)

    # create lists of means and final list of good peaks for each method and foot
    fzR_top10_data = [fzR_top10_mean, fzR_top10]
    fzR_BWcut_data = [fzR_BWcut_mean, fzR_BWcut]
    fzL_top10_data = [fzL_top10_mean, fzL_top10]
    fzL_BWcut_data = [fzL_BWcut_mean, fzL_BWcut]

    # Create nested dictionary for methods, R/L foot, mean value, peak arrays
    FinalVals = {'Top 10 Avg Method': {'Right': fzR_top10_data, 'Left': fzL_top10_data},
                 'BW Cutoff Method': {'Right': fzR_BWcut_data, 'Left': fzL_BWcut_data}}

    return FinalVals


def bar_graphs(files, FinalVals):
    N = len(files)
    files = tuple(files)

    ind = np.arrange(N)
    width = .35

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, FinalVals['Top 10 Avg Method'], width, color='b')
    rects2 = ax.bar(ind + width, FinalVals['BW Cutoff Method'], width, color='g')
    ax.set_ylabel('Trials')
    ax.set_title('Mean Forces for Right Foot')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(files)
    ax.legend((rects1[0], rects2[0]), ('Top_10R', 'BWcutR'))


def main():
    all_csv_files = []
    for root, dirs, files in os.walk('.', topdown=False):
        extensions = '.csv'
        for f in files:
            ext = os.path.splitext(f)[-1]
            if ext == extensions:
                csv_file = os.path.join(root, f)
                all_csv_files.append(csv_file)

    print('There were', str(len(all_csv_files)), '.csv files found in directory', os.path.abspath('.'))
    for c in all_csv_files:
        print(c)

    print('\n\n===================\nBegin Data Analysis\n===================\n')

    weight_file = "./data/weights.txt"
    f = open(weight_file, "r")
    weights = {}
    for line in f:
        try:
            line_split = line.split(',')
            weights[str(line_split[0])] = float(line_split[1])
        except: pass

    cleaned = {}
    for c in all_csv_files:
        csv_name = c.split('\\')[-1]
        print('--------------------')
        print(csv_name)

        if csv_name in weights:
            wgt = weights[csv_name]
            print('Located runner weight at ' + weight_file + ': ' + str(wgt) + 'kg')
        else:
            print('Error: runner weight not located in ', weight_file)
            text = 'Enter the weight of runner for file ' + str(csv_name) + ' in kg: '
            wgt = float(input(text))

        data_max = find_max(c, wgt)
        cleaned[csv_name] = data_max
        print('\nTop 10 Avg Method')
        print('Right Mean:', data_max['Top 10 Avg Method']['Right'][0])
        print('Left Mean:', data_max['Top 10 Avg Method']['Left'][0])
        print('\nBW Cutoff Method')
        print('Right Mean:', data_max['BW Cutoff Method']['Right'][0])
        print('Left Mean:', data_max['BW Cutoff Method']['Left'][0])
        print()

    #bar_graphs(all_csv_files, FinalVals)
    print(cleaned)

if __name__ == "__main__":
    main()
