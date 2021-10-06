# 2CsmPET
Script to analyze 2-color single-molecule PET intensity time traces.

The script has been developed in the programming language Python (Version 3.8) and was tested and run on a Windows operating system (Version Windows10).

Detailed instructions and explanations including generated oupput files can be found as comments in the script.

Script for analysis of two-colour single-molecule data that shows discrete fluorescence intensity
changes caused by PET fluorescence quenching. The data are analysed
in terms of dwell-times and concomitance of transitions in both
channels. A data-specific threshold is calculated
and subsequently applied to find PET-mediated fluorescence transitions beyond noise.
Traces in which transitions (steps) are found are displayed to the user who can
confirm the events, discard the trace and/or adjust default values, e.g. gaussian smoothing, for
that trace.

The script implements a code from the following link:
https://github.com/thomasbkahn/step-detect/blob/master/step_detect.py 

The script runs on python version 3.8 with specified extensions, which can be installed with e.g. pip.
Typical installation time: few minutes, depending on the installed python installation.
During execution (for example per drag and drop) the user will be asked for text-files for the RED and GREEN channels. Output
files are in the .csv file format.
For more detail information, see comments in the 2CsmPET.py file.
The python distribution Anaconda (Spyder application) can also be used to run the code. The script specifies the line of
code that must be disabled when using Spyder.
Code was executed on a Windows operating system (Version Windows10).
