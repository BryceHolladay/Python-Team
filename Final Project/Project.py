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

    plot(fzR, fzL, file.split('\\')[-1])

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


def bar_graphs(vals):
    keys = vals.keys()
    right = []
    left = []
    for k in keys:
        #r_avg = (vals[k]['Top 10 Avg Method']['Right'][0] + vals[k]['BW Cutoff Method']['Right'][0]) / 2
        #l_avg = (vals[k]['Top 10 Avg Method']['Left'][0] + vals[k]['BW Cutoff Method']['Left'][0]) / 2
        r_avg = vals[k]['Top 10 Avg Method']['Right'][0]
        l_avg = vals[k]['Top 10 Avg Method']['Left'][0]
        right.append(r_avg)
        left.append(l_avg)

    barwidth = 0.25
    data = [left, right]
    X = np.arange(len(keys))
    fig, ax = plt.subplots(figsize=(12,8))
    labels = []
    for n in range(len(data[0])):
        labels.append(data[0][n])
    for m in range(len(data[1])):
        labels.append(data[1][m])
    rects = ax.patches
    for rect, label in zip(rects, labels):
        plt.text(i.get_width()+0.2, i.get_y()+0.5,
                 labels[i])
    plt.bar(X + 0.00, data[0], color = 'b', width = barwidth, label='left')
    plt.bar(X + 0.25, data[1], color = 'r', width = barwidth, label='right')
    plt.xlabel("Trial")
    plt.ylabel('Reaction Force (N)')
    plt.xticks([r + barwidth for r in range(len(data[0]))], keys)
    for rect, label in zip(rects, labels):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2, height + 5,
                label.round(2), ha='center', va='bottom', fontweight='extra bold')

    plt.legend(loc='lower right')
    plt.title('Left and Right Force Plate Measurements')
    plt.savefig('graphs/bar_graph.png')
    print('Bar graph saved at .\\graphs\\bar_graph.png')
    return

def fft(vals):
    return

def plot(fzR, fzL, file_name):
    time = np.zeros(len(fzR)).tolist()
    t = 0
    for i in range(len(time)):
        time[i] = t
        t += 1/960
    fig = plt.figure()
    plt.plot(time, fzR, label='Right')
    plt.plot(time, fzL, label='Left')
    plt.xlabel('Time (s)')
    plt.ylabel('Reaction Force (N)')
    plt.legend(loc='upper right')
    name = 'graphs/' + file_name + '_plot.png'
    plt.savefig(name)


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
        print('Plot saved at .\\graphs\\'+ csv_name +'_plot.png\n')

    print('====================\n')
    bar_graphs(cleaned)

if __name__ == "__main__":
    main()
