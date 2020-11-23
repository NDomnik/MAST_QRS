# MAST_QRS
Automated detection of ECG R-waves, optimized for murine ECG obtained via implantable radiotelemetry

(c) Copyright 2013 Sami Torbey

DOI: 10.5281/zenodo.4287670 
__________________________________________________________________________________________________
mast_qrs is a MATLAB function that detects and returns QRS positions for a murine ECG vector input. 

Below is a MATLAB code sample that reads a text file ('ecg.txt') containing a murine ECG sampled at 100 Hz (each line representing the amplitude of one sample or 10 ms) and plots the signal, along with the detected QRS locations:


ecg = csvread('ecg.txt');

freq = 100;

qrs = mast_qrs(ecg, freq);

plot(1:length(ecg), ecg, qrs, max(ecg), 'rx')
