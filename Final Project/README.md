# ME396 Final Project
## Team Keith Lane

This project will allow you to convert raw force plate data into usable filtered data and graphs.
Data was collected by Sarah Hildreth in the UT Austin Neuromuscular Biomechanics lab group under Dr. Richard Neptune.
Participants wearing special insoles were asked to run on two treadmills, one for each foot, and Ground Reaction Forces (GRFs) were measured.
During data collection, errors were made by the participants where they placed the wrong foot on the wrong force plate during running, resulting in some data peaks reporting incorrect data.
This project exists to remove these false steps from the dataset, and return an average value of the accurately recorded peaks, as well as several modes of data visualization.
<br><br>
The Project.py script, when run, reads csv files from the data sub-directory and outputs filtered peaks and means to results.txt.
Graphs and visualizations will be placed into the graphs sub-directory.
Peaks are filtered using two methods, the Top Ten Average method and the Body Weight Cutoff method, both of which are returned and compared
Information about the process will be printed to the console for the user to review as well while the script is running.
Some plots include full graphs of the data as well as zoomed portions to better visualize the data.
<br><br>
Place new datasets into the data sub-directory, and update the file weights.txt with the user weight for any added datasets.
If no weight for a dataset is found, a weight input will be prompted from the user.
<br><br>
The graphs sub-directory contains plots of the raw data for each trial, fft results for each trial, and two bar graphs for each foot comparing the trial data to each other.
After running the script any old plots in the sub-directory will be updated, as well as the result.txt file.
