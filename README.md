# MAST_QRS
Automated detection of ECG R-waves, optimized for murine ECG obtained via implantable radiotelemetry

(c) Copyright 2011 Sami Torbey

DOI: 10.5281/zenodo.4287696 
__________________________________________________________________________________________________
mast_qrs is a MATLAB or GNU Octave function that detects QRS positions in a murine ECG vector input. If you don't already have MATLAB or GNU Octave, start by downloading and installing Octave, which is available for free on the [GNU website](https://www.gnu.org/software/octave/).

To start using MAST, we import a murine ECG time series into a GNU Octave or MATLAB variable. Multi-lead ECGs can be represented as matrices where rows correspond to time and each column corresponds to amplitude for one lead. In the example below, we read the time series from a text file ('ecg.txt') containing one sample per line:

    ecg = dlmread('ecg.txt');

Next we exclude any leads we don't want used for QRS detection. Here our 'ecg' variable contains multiple leads and we only want to use leads 2 and 3 for QRS detection:

    ecg = ecg(:, [2, 3]);

Now we're ready to run MAST:

	[qrs, valid] = mast_qrs(ecg, 1000);

This detects QRS locations in the 'ecg' time series and saves the indices of QRS midpoints in the 'qrs' variable. The second parameter (1000 in this example) is optional and indicates the sampling frequency in Hz (if omitted, mast_qrs will assume a 1,000 Hz sampling frequency).

MAST automatically excludes detections that it considers highly likely to be false positives. The second return value ('valid' in this example) is optional, and indicates the range of indices that MAST considered valid (hence not excluded). We recommend limiting any analysis to these valid ranges.

Lastly, let's plot the input ECG signal along with MAST's QRS detections (red vertical lines) and the valid ranges (green horizontal lines):

    plot(1:length(ecg), ecg, '', repmat(qrs,1,2), [min(ecg(:)) max(ecg(:))], 'r', valid', repmat(max(ecg(:)),1,2), 'g')
