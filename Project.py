## ME396P Keith Lane Final Project
## Team Members: Ali Bhaiwala, Sarah Hildreth, Bryce Holladay
## Spring 2021

from scipy.signal import find_peaks
from scipy import fft
import pandas as pd
import matplotlib.pyplot as plt

# Function that takes in a csv file and outputs an average max peak, as well as saves a graph
def find_max(file):
    csv = pd.read_csv(file, header=[3], nrows=36000) #Header [3] to only include 4th row of column titles
        
    fzR = csv['Fz'].tolist() #creating a list of all Fz values on right foot
    fzL = csv['Fz.1'].tolist() #create list of all Fz values on left leg
       
    return

find_max('S4_T21.csv')
