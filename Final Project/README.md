# ME396 Final Project
## Team Keith Lane

This directory will allow you to convert raw force plate data into usable filtered data and graphs.
The Project.py script, when run, reads csv files from the data sub-directory, and outputs data to results.txt and to the graphs sub-directory.
Information about the process will be printed to the console for the user to review as well.
Some plots include full graphs of the data as well as zoomed portions to better visualize the data.
<br><br>
Place new datasets into the data sub-directory, and update the file weights.txt with the user weight for this dataset.
If no weight for a dataset is found, a weight will be prompted from the user.
<br><br>
The graphs sub-directory contains plots of the raw data for each trial, fft results for each trial, and two bar graphs for each foot comparing the trial data to each other.
After running the script any old plots in the sub-directory will be updated.