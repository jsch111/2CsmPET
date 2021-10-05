# 2CsmPET
Script to analyze 2-color single-molecule PET intensity time traces.

Script for analysis of two-color single-molecule data which shows digital intensity
changes due to PET fluorescence fluctuations. The data are analyzed
in terms of dwell-times and concomittance of transitions between both
channels, thereby iterating over traces. A data-specific threshold is calculated
and subsequently applied to find PET-mediated transitions and to ignore noise.
Traces in which transitions (steps) are found are dispalyed to the user who can
confirm the events, discard the trace and/or change default values for
that trace.

The script uses code from the following link:
https://github.com/thomasbkahn/step-detect/blob/master/step_detect.py 

The script runs on a working version of python 3.8 with the specified extensions, which can be installed with pip, for example.
During execution (for example per drag and drop) the user will be asked for text-files for the RED and GREEN channels. 
For more information, see the 2CsmPET.py file.
The python distribution Anaconda (Spyder application) can also be used to run the code. The script specifies the line of
code that must be disabled when using Spyder.
