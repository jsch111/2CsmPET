"""
Script to analyze 2-color single-molecule PET intensity time traces.

Usage
----------
Script for analysis of two-color single-molecule data which shows digital intensity
changes due to PET fluorescence fluctuations. The data are analyzed
in terms of dwell-times and concomittance of transitions between both
channels, thereby iterating over traces. A data-specific threshold is calculated
and subsequently applied to find PET-mediated transitions and to ignore noise.
Traces in which transitions (steps) are found are dispalyed to the user who can
confirm the events, discard the trace and/or change default values for
that trace.

Notes
----------
Data-specifc threshold:
A default threshold of 0.01 is applied. Function iterates over traces and several
normalization and resizing steps of the data and substraction of mean intensity is 
followd by smoothing. The multiscale product is calculated using the step_detect.py
module from Thomas Kahn. Further, the module is used to find steps above the initial
threshold. Relative magnitudes ("peak-heights") of all steps are collected with 
the step_detect module. The average value is calculated and serves subsequently as
data-specific threshold.

definition (bleaching):
Frame counts as a belaching event if the signal remains eblow the average value
from frame 5 to the current frame for longer than the number of frames specified
by belaching_index. Set by default to 341, representing 340 frames.

concomittance:
A transition is identified as simultaneous if the events occurs in both channels
within the indicated time window (default: 6 frames).

calculation of events:
Traces (events) are filtered to exhibit any intensity fluctuations in both channels
to filter for functional molecules. Further, only alternating transitions are
considered as these originate from PET. Events are sorted in simultaneous and
not simultaneous ones by considering the concomittance time-window.

delete events option:
During analysis events can be manually selected to be not considered. Reasons for this
are: transition with low step magnitude (being not caused by PET) incorrectly identified 
as step due to insufficient thresholding; Steps that occur after a bleaching event,
that was not detected properly; Incorrectly detected ON-OFF alternation (instead of
OFF-OFF or ON-ON) due to non-detection of a step beforehand.

supported Files:
text Files, with spacebar as delimiter

output-File:
CSV-File, including ON and OFF dwell times for each channel as well as the absolute
number of transitions/steps that occured simultaneous and not simultaneous.


References
----------
Thomas Kahn; step_detect.py module from:
https://github.com/thomasbkahn/step-detect/blob/master/step_detect.py 

Brian M. Sadler, Ananthram Swami (1998): Analysis of wavelet transform multiscale products
for step detection and estimation. In: Army Research Laboratory Tech. Rep. ARL-TR-1664

Sadler, B. M.; Swami, A. (1999): Analysis of multiscale products fo step detection and
estimation. In: IEEE Trans. Inform. Theory 45 (10), S. 630-642. DOI: 10.1038/nrm3658

Author: Jonathan Schubert
Licence: MIT
"""

import os
import time
from itertools import zip_longest
from copy import deepcopy
import csv
import warnings

from tkinter import filedialog
from tkinter import Tk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import step_detect
from scipy.ndimage import gaussian_filter1d


warnings.simplefilter(action = "ignore", category = RuntimeWarning)

root = Tk()
root.withdraw()


def default_thresholds():
    """	Default values.
	
	sigma_r/g = int
	            standard deviation for gaussian kernel of 
                gaussian_filter1d for red and green channel, set 5 by default.
	threshold_r/g = float
	                initial threshold used for calculation of data specific threshold, set 0.01 by default.
	n = int
	    scales to multiply to, used by mz_fwt. Set to 3 by default. If n is odd directionality of the step
        is preserved, if n is even directionality gets lost. For more information refer to step_detect.py.
	delete_red, delete_green = list
	                           list of events, that should be not considered for that trace due to origins other than PET.
		                       for more details refer to "Notes".
    
	Returns: Default values.
    """
    sigma_r = 5
    sigma_g = 5
    threshold_r = 0.01
    threshold_g = 0.01
    n = 3
    delete_red = "Type the red events to be deleted (first event = 1, ...)"
    detele_green = "Type the green events to be deleted (first event = 1, ...)"
    return sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, detele_green


    
def change_thresholds(sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, delete_green):
    """	Function to change thresholds used in default_thresholds() while script is running.
    In addition the option to not consider events manually is provided. For more details
	refer to "Notes".
    
    Parameters are loaded from default_thresholds().
    
    Returns updated parameters.
    """
    options = [sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, delete_green]
    dic_int = ["sigma_red", "sigma_green", "threshold_red", "threshold_green", "n (scale to multiply to)", "delete_red", "delete_green"]
    while True:
        try:
            th = int(input("Which threshold do you want to change?\nsigma_red = 1\nsigma_green = 2\nthreshold_red = 3\nthreshold_green = 4\nn (scale to multiply to) = 5\ndelete_red = 6\ndelete_green = 7\n"))
            break
        except:
            print("Try again and use only numbers between 1 and 7!")
    if th < 6:
        for i in range(5):
            if i == th-1:
                while True:
                    try:
                        options[i] = float(input("(Old = %s) New %s = " %(options[i], dic_int[i])))
                        print("New threshold: %s = %s" %(dic_int[i], options[i]))
                        break
                    except:
                        print("Try again and please use only numbers!")
    elif th == 6:
        while True:
            try:
                de = input("Type the red event(s) to be deleted (separated by a comma, e.g. '1, 4'): ")
                dele = de.split(",")
                options[5] = [int(x) for x in dele]
                print(options[5])
                break
            except:
                print("Try again and please use only one or more NUMBERS separated by a comma!")
    elif th == 7:
        while True:
            try:
                de = input("Type the green event(s) to be deleted (separated by a comma, e.g. '1, 4'): ")
                dele = de.split(",")
                options[6] = [int(x) for x in dele]
                print(options[6])
                break
            except:
                print("Try again and please use only one or more NUMBERS separated by a comma!")
    else:
        print("Try again and use only numbers betwenn 1 and 7!")
    
    sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, delete_green = options[0], options[1], options[2], options[3], options[4], options[5], options[6]
    return sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, delete_green
     


def delete_event_red(steps, delete_red):
    """	Performs ignoration of incorrectly detected events if the option 
    in change_thresholds() is selected and events were typed in for red channel.
    
    Parameters:
    -----------
    steps = list
            List of indices of steps.
    
    Returns:
    --------
    steps = list
            Updated List of indices.
    """
    if not type(delete_red) == str:
        delete_red.reverse()
        for i in delete_red:
            del steps[i-1]
    return steps



def delete_event_green(steps, delete_green):
    """	
    Performs ignoration of incorrectly detected events if the option 
    in change_thresholds() is selected and events were typed in for green channel.
    
    Parameters:
    -----------
    steps = list
            List of indices of steps.
    
    Returns:
    --------
    steps = list
            Updated List of indices.
    """
    if not type(delete_green) == str:
        delete_green.reverse()
        for i in delete_green:
            del steps[i-1]
    return steps



def step_sizes_own(data, indices):
    """ Determines the directionality of each step (index) by looking up
    the values of the input data (mz_fwt from step_detect.py) and sorting them
    into negative (<0) and positive (>0).
    
    Parameter:
    ----------
    data = numpy array
           The multiscale product of data
    indices = list
              List of indices of the detected steps (see step_detect.py)
    
    Returns:
    --------
    peaks_pos, peaks_neg = list
                           Lists of indices of positive and negative steps, respectively.
    """
    peaks_pos = []
    peaks_neg = []
    indices = sorted(indices)
    for i in indices:
        if data[i] > 0:
            peaks_pos.append(i)
        elif data[i] < 0:
            peaks_neg.append(i)
    peaks_neg = [i*-1 for i in peaks_neg]
    return peaks_pos, peaks_neg


def calculate_threshold(raw_data, sigma, threshold):
    """ Calculation of data-specific threshold.
    
    Parameters:
    -----------
    raw_data = list
               Contains raw data set after loading it into the script. Time series of data points.
    sigma = scalar
            Standard deviation for gaussian kernel. (for further 
            information refer to gaussian_filter1d, numpy).
    threshold = float
                Threshold value that has to be reached by peaks to be considered.

    Returns:
    --------
    calculated_threshold = float
                           Data specific threshold calculated by averaging of relative magnitudes
                           of all peaks/steps that reached initial threshold.
    """
    step_sizes = []
    for i in range(raw_data.shape[1]-1):    
        sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, delete_green = default_thresholds()
        i = i+1                             
        data = np.array(raw_data[:, i])
        data /= np.abs(data).max()
        data -= np.average(data)
        data_smoothed = gaussian_filter1d(data, sigma, order=0)
        data_smoothed /= np.abs(data_smoothed).max()
        mz_fwt = step_detect.mz_fwt(data_smoothed, n=n)
        mz_fwt *= 100
        mz_fwt /= np.abs(mz_fwt).max()
        mz_fwt = np.roll(mz_fwt, sigma-2)
        steps = step_detect.find_steps(mz_fwt, threshold)
        peaks_pos, peaks_neg = step_sizes_own(mz_fwt, steps)
        step_size = []
        step_size = [mz_fwt[i] for i in steps]
        step_sizes.extend(np.abs(step_size))
    calculated_threshold = np.mean(step_sizes)
    print("calculated threshold :", calculated_threshold)
    return calculated_threshold




def check_for_bleaching(red, green, bleaching_index=341):
    """ Searches for bleaching events in smoothed data. 
    
    Parameters:
    -----------
    red/green = numpy array
                1-dimensional array that represents time series of datapoints
    bleaching_index = int
                      Number of frames, for which the signal has to remain below the average 
					  value from frame 5 to the current one, to be detectes as bleaching event.
    
    Returns:
    --------
    list1/list2 = list 
                  Lists contain frame of the bleaching event in the respective data channel.
    """
	#red
    list1 = []
    for idx, i in enumerate(red):   
       a = idx
       b = red[idx:idx+bleaching_index]
       if np.all(b < np.mean(red[5:idx])) == True:
           if not list1:                # checks if list1 is empty, if yes it proceeds
               print("RED bleaching event at frame", a)
               list1.append(a)
           else:
               pass
       else:
           pass
    #green
    list2 = []
    for idx2, i2 in enumerate(green):   
       a2 = idx2
       b2 = green[idx2:idx2+bleaching_index]
       if np.all(b2 < np.mean(green[5:idx])) == True:
           if not list2:                # checks if list2 is empty, if yes it proceeds
               print("GREEN bleaching event at frame", a2)
               list2.append(a2)
           else:
               pass
       else:
           pass
       
    return list1, list2



def consider_bleaching(bleaching_event_red, bleaching_event_green, times_copy, times_g_copy):
    """ Considers present bleaching events. If bleaching events were found by check_for_bleaching()
    following events are not considered. Returns updated lists of steps.
	
	Parameters:
	-----------
	bleaching_event_red, bleaching_event_green = list
                                                 contains the frame of the bleaching event
	times_copy, times_g_copy = list
	                           Copy of lists of indices of steps

	Returns:
	--------
	times, times_g = list
	                 Updated version of lists of indiced of steps
    """
    times
    times_g
    if bleaching_event_red:
        for value in times_copy:
            if abs(value) >= bleaching_event_red[0]:
                times.remove(value)
            else:
                pass
        for value_g in times_g_copy:
            if abs(value_g) >= bleaching_event_red[0]:
                times_g.remove(value_g)
            else:
                pass
    else:
        pass
    
    if bleaching_event_green:
        for value_g2 in times_g_copy:
            if abs(value_g2) >= bleaching_event_green[0]:
                times_g.remove(value_g2)
            else:
                pass
        for value2 in times_copy:
            if abs(value2) >= bleaching_event_green[0]:
                times.remove(value2)
            else:
                pass
    else:
        pass
    return times, times_g




def simultancy(red, green, window):
    """ Checks for concomitance of steps/events in red and green channel and divide
    them in simultaneous and not simultaneous ones. Events are sorted for absolute value
    and gets compared iteratively against each other to detect steps that occur 
    simultaneous within the selected time window.
    
    Parameters:
    -----------
    red/green = list
                List of indices of steps.
    window = int
             Number of frames that is used as time window to divide events
             into simultaneous and not simultaneous ones. Starting from 
             red event, code is searching for event in green channel that occurs within 
             time window before or after the frame/time point of interest and vice versa.
             E.g. window = 5, Event in red channel at frame = 51. Window for green 
             (and red) event to be assigned as simultaneous = frame 47-55.
    
    Returns:
    --------
    simul_r/simul_g = list
                      Lists contain events/steps that occur simultaneous in the particular channel.
    non_simul_r/non_simul_g = list
                              Lists contain events/steps that occur not simultaneous in the particular channel.
    """
    window = window-1       # actual window size equal to 'window'
    red_events = red
    red_events.sort(key=abs)
    green_events = green
    green_events.sort(key=abs)

    simul_r = []
    non_simul_r = []
    simul_g = []
    non_simul_g = []

    counter = 1
    while counter < len(red_events):
        for event_r in red_events:
            list_g = []
            for event_g in green_events:
                for z in range(event_g, event_g+window+1):
                    list_g.append(z)
                for z2 in range(event_g-window, event_g):
                    list_g.append(z2)
            if event_r in list_g:
                simul_r.append(event_r)
                counter += 1
            else:
                non_simul_r.append(event_r)
                counter += 1
    print("red simul/not simul: ", simul_r, "/", non_simul_r)
    
    counter2 = 1
    while counter2 < len(green_events):
        for event_g in green_events:
            list_r = []
            for event_r in red_events:
                for y in range(event_r, event_r+window+1):
                    list_r.append(y)
                for y2 in range(event_r-window, event_r):
                    list_r.append(y2)
            if event_g in list_r:
                simul_g.append(event_g)
                counter2 += 1
            else:
                non_simul_g.append(event_g)
                counter2 += 1
    print("green simul/not simul: ", simul_g, "/", non_simul_g)

    return simul_r, non_simul_r, simul_g, non_simul_g




def calculate_events(red, green):
    """ Calculates step-durations (dwell times) for both channels using indices of detected steps. 
    Filter is implemented that allows only events that occur alternatingly (ON-OFF or OFF-ON). 
    In order to analyze only functional molecules only traces (their steps in particular) 
    are allowed that show minimum one ON- and OFF-dwell time per channel. Otherwise trace is skipped. 
	
	Parameters:
	-----------
	red, green = list
	             Lists of indices of detected steps of the respective channel after considered bleaching events
	
	Returns:
	--------
	on_peaks_red, off_peaks_red = list
	                              contains ON or OFF steps found in red channel
	on_peaks_green, off_peaks_green = list
	                                  contains ON or OFF steps found in green channel
	results_on1, results_off1 = list
	                            contains dwell times of ON or OFF states, red channel
	results_on2, results_off2 = list
	                            contains dwell times of ON or OFF states, green channel
	times, times_g = list
	                 contains unfiltered steps (after consideres bleaching events) of red or green channel
	peaks, peaks_g = list
	                 updated version of times, times_g
    """
    times = red
    times_g = green
    
    #red; checks for alternating ON/OFF events + calculation of dwell time
    times.sort(key=abs)
    results_on1 = []
    results_off1 = []
    for j in range(len(times)-1):
        result1 = abs(times[j+1]) - abs(times[j])
        if times[j+1] > 0 and times[j] > 0:
            pass
        elif times[j+1] < 0 and times[j] < 0:
            pass
        elif times[j+1] > 0:
            results_off1.append(result1)
        elif times[j+1] < 0:
            results_on1.append(result1)
    
    #green; checks for alternating ON/OFF events + calculation of dwell time
    times_g.sort(key=abs)
    results_on2 = []
    results_off2 = []
    for j2 in range(len(times_g)-1):
        result2 = abs(times_g[j2+1]) - abs(times_g[j2])
        if times_g[j2+1] > 0 and times_g[j2] > 0:
            pass
        elif times_g[j2+1] < 0 and times_g[j2] < 0:
            pass
        elif times_g[j2+1] > 0:
            results_off2.append(result2)
        elif times_g[j2+1] < 0:
            results_on2.append(result2)
    
    # results are only considered if minimum 1 ON- and OFF-dwell time per channel present
    if results_on1 and results_off1 and results_on2 and results_off2:
        #red
        on_peaks_red = [i for i in times if i>0]
        off_peaks_red = [i for i in times if i<0]
        peaks = on_peaks_red + off_peaks_red
        peaks.sort(key=abs)
        #green
        on_peaks_green = [i for i in times_g if i>0]
        off_peaks_green = [i for i in times_g if i<0]
        peaks_g = on_peaks_green + off_peaks_green
        peaks_g.sort(key=abs)
    else:
        # if less than one event per channel lists are cleared
        results_on1 = []
        results_off1 = []
        results_on2 = []
        results_off2 = []
    
    return on_peaks_red, off_peaks_red, results_on1, results_off1, on_peaks_green, off_peaks_green, results_on2, results_off2, peaks, peaks_g, times, times_g


 
def plot_data(data_o_r, data_o_g, peaks_red, peaks_green, data_smoothed_r, data_smoothed_g, sigma, orig_data_red, orig_data_green):
    """ Function to plot different data analyzing steps in order to confirm results manually while script is running.
    4 horizontally ordered panels. Found steps are superimposed as dotted vertical lines in each panel. In panel
    1 and 2 red and green colored bars indicate ON states (light colors) and OFF states (dark colors).
        Panel 1: raw intensity time traces of red and green channel are plotted
        Panel 2: data after substraction of average intensities and normalization to 1
        Panel 3: data smoothed with gaussian_filter1d
        Panel 4: multiscale product of mz discrete forward wavelet transform, calculated with mz_fwt from step_detect.py.
                 Calculated thresholds for red and green data are indicated as horizontal colored lines.

    Parameters:
    -----------
    data_o_r, data_o_g = numpy array
                         1-dimensional array that represents time series of datapoints, red or green channel
						 After smoothing and substraction of average
    data_smoothed_r, data_smoothed_g = numpy array
                                       1-dimensional array that represents time series of datapoints, red and green channel
                                       After smoothing with gaussian_filter1d and normalization to 1
    orig_data_red, orig_data_green = numpy array
                                     1- dimensional array that represents time series of datapoints, red or green channel
									 raw intensity time traces
	peaks_red, peaks_green = list
	                         contain steps that will be considered, red or green channel
    """
    threshold_r
    threshold_g
    peaks_red_abs = []
    peaks_green_abs = []
    
    data_o_r /= np.abs(data_o_r).max()   # for better view
    data_o_g /= np.abs(data_o_g).max()   # for better view

    for value in peaks_red:
        peaks_red_abs.append(abs(value))
    for value2 in peaks_green:
        peaks_green_abs.append(abs(value2))
    style.use("seaborn-whitegrid")
    f,plt_arr = plt.subplots(4, sharex=True)
    f.suptitle("trace %i" %(i))
    plt_arr[0].plot(orig_data_red, color="red", alpha=0.7)
    plt_arr[0].plot(orig_data_green, color="green", alpha=0.7)
    plt_arr[0].set_title("raw intensity trace")
    plt_arr[0].set_ylabel("fluorescence intensity", fontsize="small")
    for ii in peaks_red_abs:
        plt_arr[0].axvline(x=ii, linestyle=":", color="red")
    for jj in peaks_green_abs:
        plt_arr[0].axvline(x=jj, linestyle="--", color="green", alpha=0.7)
        
    for kk in range(len(peaks_red)-1):
        if peaks_red[kk] < 0 and peaks_red[kk+1] > 0:
            plt_arr[0].hlines(y=(orig_data_red.max()), xmin=abs((peaks_red[kk])), xmax=abs((peaks_red[kk+1])), linewidth=5, color="darkred")
        elif peaks_red[kk] > 0 and peaks_red[kk+1] < 0:
            plt_arr[0].hlines(y=(orig_data_red.max()), xmin=abs((peaks_red[kk])), xmax=abs((peaks_red[kk+1])), linewidth=5, color="red")
        else:
            pass
        
    for ll in range(len(peaks_green)-1):
        if peaks_green[ll] < 0 and peaks_green[ll+1] > 0:
            plt_arr[0].hlines(y=(orig_data_red.max()-10), xmin=abs((peaks_green[ll])), xmax=abs((peaks_green[ll+1])), linewidth=5, color="darkgreen")
        elif peaks_green[ll] > 0 and peaks_green[ll+1] < 0:
            plt_arr[0].hlines(y=(orig_data_red.max()-10), xmin=abs((peaks_green[ll])), xmax=abs((peaks_green[ll+1])), linewidth=5, color="lime")
        else:
            pass

    plt_arr[0].set_ylabel("fluorescence intensity", fontsize="small")
    plt_arr[0].set_xlim(0, 1800)

    plt_arr[1].plot(data_o_r, color="red", alpha = 0.5)
    plt_arr[1].set_title("normalized fluorescence intensity")
    plt_arr[1].plot(data_o_g, color="green", alpha = 0.5)
    for ii in peaks_red_abs:
        plt_arr[1].axvline(x=ii, linestyle=":", color="red")
    for jj in peaks_green_abs:
        plt_arr[1].axvline(x=jj, linestyle="--", color="green", alpha=0.7)

    for kk in range(len(peaks_red)-1):
        if peaks_red[kk] < 0 and peaks_red[kk+1] > 0:
            plt_arr[1].hlines(y=0.90, xmin=abs((peaks_red[kk])), xmax=abs((peaks_red[kk+1])), linewidth=5, color="darkred")
        elif peaks_red[kk] > 0 and peaks_red[kk+1] < 0:
            plt_arr[1].hlines(y=0.90, xmin=abs((peaks_red[kk])), xmax=abs((peaks_red[kk+1])), linewidth=5, color="red")
        else:
            pass
        
    for ll in range(len(peaks_green)-1):
        if peaks_green[ll] < 0 and peaks_green[ll+1] > 0:
            plt_arr[1].hlines(y=0.85, xmin=abs((peaks_green[ll])), xmax=abs((peaks_green[ll+1])), linewidth=5, color="darkgreen")
        elif peaks_green[ll] > 0 and peaks_green[ll+1] < 0:
            plt_arr[1].hlines(y=0.85, xmin=abs((peaks_green[ll])), xmax=abs((peaks_green[ll+1])), linewidth=5, color="lime")
        else:
            pass
        
    plt_arr[1].set_ylim(-1, 1)


    
    plt_arr[2].plot(data_smoothed_r, color="red", alpha = 0.5)
    plt_arr[2].plot(data_smoothed_g, color="green", alpha = 0.5)
    plt_arr[2].set_title("smoothed with gaussian_filter1d, sigma = %i" %(sigma))
    for ii in peaks_red_abs:
        plt_arr[2].axvline(x=ii, linestyle=":", color="red")
    for jj in peaks_green_abs:
        plt_arr[2].axvline(x=jj, linestyle="--", color="green", alpha=0.7)

    plt_arr[3].plot(mz_fwt, color="red", alpha = 0.5)
    plt_arr[3].plot(mz_fwt_g, color="green", alpha = 0.5)
    plt_arr[3].set_xlabel("frames")
    plt_arr[3].set_title("Multiscale product of the MZ discrete forward wavelet transform")
    plt_arr[3].axhline(y=threshold_r, linestyle="solid", linewidth=1.0, color="red")
    plt_arr[3].axhline(y=(threshold_r*-1.0), linestyle="solid", linewidth=1.0, color="red")
    plt_arr[3].axhline(y=threshold_g, linestyle="solid", linewidth=1.0, color="lime")
    plt_arr[3].axhline(y=(threshold_g*-1.0), linestyle="solid", linewidth=1.0, color="lime")
    for ii in peaks_red_abs:
        plt_arr[3].axvline(x=ii, linestyle=":", color="red")
    for jj in peaks_green_abs:
        plt_arr[3].axvline(x=jj, linestyle="--", color="green", alpha=0.7)
    
    mng = plt.get_current_fig_manager()         # maximizes actual window
    mng.window.state("zoomed")                  # maximizes actual window; has to be deactivated, if spyder is used to run the code
    plt.show(block=False)
    


	
""" Script asks for two .txt files that represent RED and GREEN time series and files are loaded
using numpy. Files get reduced to every second column (contain maximum grey values).
After that different functions defined above are used to analyze the red and green data
sets. Therefor script iterates over traces in red and green channel simultaneous.
Also steps that are detailly described in the calculated_threshold function are
applied to the data thereby using the calculated threshold.
Traces in which steps were found are shown to the user with the plot_data function
that allows user to change initial values for the particular trace (see default_thresholds)
and confirmation of the found steps/dwell times.
"""

File_red = filedialog.askopenfilename(initialdir = "D:/Promotion/1_single-molecule PET fluorescence imaging microscopy/2-Farben-Messungen",
                                      title = "Select input file for RED dye ", filetypes = (("text files", "*.txt"), ("all files", "*.*")))
FileName = os.path.splitext(os.path.basename(File_red))[0]
DirName = os.path.dirname(File_red)

File_green = filedialog.askopenfilename(initialdir = DirName,
                                        title = "Select input file for GREEN dye ", filetypes = (("text files", "*.txt"), ("all files", "*.*")))
FileName = os.path.splitext(os.path.basename(File_green))[0]
DirName = os.path.dirname(File_green)

start = time.time()

content_r = np.loadtxt(File_red)
content_g = np.loadtxt(File_green)
col_redu = len(content_g)
red_raw = content_r[:col_redu, ::2]
green_raw = content_g[:, ::2]


sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, delete_green = default_thresholds()
print("RED: ")
calculated_threshold_red = calculate_threshold(red_raw, sigma_r, threshold_r)
print("GREEN: ")
calculated_threshold_green = calculate_threshold(green_raw, sigma_g, threshold_g)


red_on_steps = []
red_off_steps = []
red_on_durations = []
red_off_durations = []
green_on_steps = []
green_off_steps = []
green_on_durations = []
green_off_durations = []

simultaneous_events_red = []
simultaneous_events_green = []
non_simultaneous_events_red = []
non_simultaneous_events_green = []

for (i, k) in zip(range(red_raw.shape[1]-1), range(green_raw.shape[1]-1)):     
    sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, delete_green = default_thresholds()
    threshold_r = calculated_threshold_red
    threshold_g = calculated_threshold_green
    i = i+1                         # skips first column (numeration)
    k = k+1                         # skips first column (numeration)
    v = 0
    while v != 1:
        try:
            print("-------------\nTrace %s" %(i))
            orig_data_red = np.array(red_raw[:, i])        # for plot_data option of raw traces
            data_o_r = np.array(red_raw[:, i])
            data_o_r /= np.abs(data_o_r).max()
            data_o_r -= np.average(data_o_r)
            data_smoothed_r = gaussian_filter1d(data_o_r, sigma_r, order=0)
            data_smoothed_r /= np.abs(data_smoothed_r).max()
            mz_fwt = step_detect.mz_fwt(data_smoothed_r, n=n)
            mz_fwt *= 100
            mz_fwt /= np.abs(mz_fwt).max()
            mz_fwt = np.roll(mz_fwt, sigma_r-2)
            steps = step_detect.find_steps(mz_fwt, threshold_r)
            steps = delete_event_red(steps, delete_red)
            peaks_pos, peaks_neg = step_sizes_own(mz_fwt, steps)
            
            # green channel
            orig_data_green = np.array(green_raw[:, k])        # for plot_data option of raw traces
            data_o_g = np.array(green_raw[:, k])
            data_o_g /= np.abs(data_o_g).max()
            data_o_g -= np.average(data_o_g)
            data_smoothed_g = gaussian_filter1d(data_o_g, sigma_g, order=0)
            data_smoothed_g /= np.abs(data_smoothed_g).max()
            mz_fwt_g = step_detect.mz_fwt(data_smoothed_g, n=n)
            mz_fwt_g *= 100
            mz_fwt_g /= np.abs(mz_fwt_g).max()
            mz_fwt_g = np.roll(mz_fwt_g, sigma_g-2)
            steps_g = step_detect.find_steps(mz_fwt_g, threshold_g)
            steps_g = delete_event_green(steps_g, delete_green)
            peaks_pos_g, peaks_neg_g = step_sizes_own(mz_fwt_g, steps_g)
            
            # bleaching-test; 1. checks for bleaching-events
            red_bleaching, green_bleaching = check_for_bleaching(data_o_r, data_o_g)
            times = peaks_pos + peaks_neg
            times.sort(key=abs)
            times_g = peaks_pos_g + peaks_neg_g
            times_g.sort(key=abs)
            times_copy = deepcopy(times)
            times_g_copy = deepcopy(times_g)
            
            # bleaching-test; 2. if bleaching event is present, events happening afterwards gets deleted
            times, times_g = consider_bleaching(red_bleaching, green_bleaching, times_copy, times_g_copy)

            on_peaks_red, off_peaks_red, results_on1, results_off1, on_peaks_green, off_peaks_green, results_on2, results_off2, peaks, peaks_g, times, times_g = calculate_events(times, times_g)
            if peaks and peaks_g:
                plot_data(data_o_r, data_o_g, peaks, peaks_g, data_smoothed_r, data_smoothed_g, sigma_r, orig_data_red, orig_data_green)
            
                print("red events to be added: ", peaks)
                print("green events to be added: ", peaks_g)
                
                choice = input("Do you agree? Press ENTER. If not press 'x' to skip this trace or 't' to modify thresholds: ")
                if choice == "x":
                    plt.close()
                    print("Trace", i, "skipped!")
                    break
                elif choice == "t":
                    sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, delete_green = change_thresholds(sigma_r, sigma_g, threshold_r, threshold_g, n, delete_red, delete_green)
                    plt.close()
                    v = 0
                else:
                    red_on_steps.extend(on_peaks_red)
                    red_off_steps.extend(off_peaks_red)
                    red_on_durations.extend(results_on1)
                    red_off_durations.extend(results_off1)
                    green_on_steps.extend(on_peaks_green)
                    green_off_steps.extend(off_peaks_green)
                    green_on_durations.extend(results_on2)
                    green_off_durations.extend(results_off2)
                              
                    simul_r, non_simul_r, simul_g, non_simul_g = simultancy(times, times_g, 6)
                    simultaneous_events_red.extend(simul_r)
                    simultaneous_events_green.extend(simul_g)
                    non_simultaneous_events_red.extend(non_simul_r)
                    non_simultaneous_events_green.extend(non_simul_g)
                    print("current steps RED", red_on_steps, red_off_steps)
                    print("current steps GREEN ", green_on_steps, green_off_steps)
                    plt.close()
                    v = 1
            else:
                print("---nothing found---", "less than 1 event per channel")
        except:
            print("ERROR: nothing found")
            plt.close()
            break
        plt.close()

red_events_simul = [len(simultaneous_events_red)]
green_events_simul = [len(simultaneous_events_green)]
red_events_not = [len(non_simultaneous_events_red)]
green_events_not = [len(non_simultaneous_events_green)]

print("-----------------------------------")
print("red_steps ON: ", red_on_steps)
print("red_steps OFF: ", red_off_steps)
print("green_steps ON: ", green_on_steps)
print("green_steps OFF: ", green_off_steps)
print("red_step ON-durations: ", red_on_durations)
print("red_step OFF-durations: ", red_off_durations)
print("green_step ON-durations: ", green_on_durations)
print("green_step OFF-durations: ", green_off_durations)
print("-----------------------------------")
print("ON counts (red):", len(red_on_durations))
print("OFF counts (red):", len(red_off_durations))
print("RED simultaneous events:", len(simultaneous_events_red))
print("RED not simultaneous events:", len(non_simultaneous_events_red))
print("ON counts (green):", len(green_on_durations))
print("OFF counts (green):", len(green_off_durations))
print("GREEN simultaneous events:", len(simultaneous_events_green))
print("GREEN not simultaneous events:", len(non_simultaneous_events_green))
print("-----------------------------------")
end = time.time()
print("--- %s seconds ---" %(end -start))
print("-----------------------------------")


out = [red_on_durations, red_off_durations, green_on_durations, green_off_durations,
       red_events_simul, red_events_not, green_events_simul, green_events_not]
outputpath = "%s\\%s_out_own.csv" %(DirName, FileName)
outputlist = list(zip_longest(*out, fillvalue=""))

with open(outputpath, "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["ON red", "OFF red", "ON green", "OFF green", "simultaneous events RED", 
                     "not simultaneous events RED", "simultaneous events GREEN", 
                     "not simultaneous events GREEN"])
    writer.writerows(zip_longest(*out, fillvalue=""))


print("End of File")
input()
    

