% MAST_QRS: Automated detection of ECG R-waves, optimized for murine ECG 
% obtained via implantable radiotelemetry
% Copyright (C) 2011 Sami Torbey

% Input:
%   ecg: single-lead or multi-lead ECG time series (leads are columns)
%   freq: ECG sampling rate (default: 1000 samples per second)
% Outputs:
%   qrs_locations: indices of detected QRS locations in ecg
%   valid_ranges: ranges of indices where QRS detection is deemed trustworthy
%                 (default: automatically drop untrustworthy ranges)
% Usage example:
%   % Return trustworthy detected QRS indices for time series "ecg"
%   qrs = mast_qrs(ecg);
%   % Plot time series "ecg" along with vertical red lines depicting detected
%   % QRS locations
%   plot(1:length(ecg),ecg,'',[qrs qrs],[min(ecg(:)) max(ecg(:))],'r')

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU LGPL        Public License as published by the Free
% Software Foundation, version 3 of the License.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU LGPL 
% Lesser General Public License for more details.

% You should have received a copy of the GNU LGPL Lesser General Public 
% License along with this program. It is available as file <COPYING.LESSER>
% in the repository alongside GNU license <LICENSE>.
% If not, or for further information, see <http://www.gnu.org/licenses/#LGPL>

function [ qrs_locations, valid_ranges ] = mast_qrs( ecg, freq )

if (nargin < 2),
    freq = 1000;
end

% Configure parameter values
std_coef = 3; % detection sensitivity parameter
narrow_window = round(0.01*freq); % window of length close to an average QRS
wide_window = round(0.1*freq); % window of length close to an average RR interval
valid_window = round(0.01*freq); % window for flagging untrustworthy detections
valid_threshold = 3/4; % threshold for flagging untrustworthy detections
merge_window = round(0.003*freq); % merge points within a single QRS

% Subtract baseline and combine leads
rectified = abs(ecg-movmean(ecg,5*wide_window));
rectified = rescale(rectified, 'InputMin', min(rectified), 'InputMax', max(rectified));
rectified = sum(rectified,2);

% Compute feature and threshold signals
narrow_threshold = movmean(rectified,narrow_window) ...
    + std_coef .* movstd(rectified,narrow_window);
wide_threshold = movmean(rectified,wide_window) ...
    + std_coef .* movstd(rectified,wide_window);

% Find candidate QRS locations and refine them
candidate_locations = find(narrow_threshold>wide_threshold);
diff_locations = find(diff([-1; candidate_locations])>1);
qrs_locations = [];
if (~isempty(diff_locations))
    qrs_locations = candidate_locations([diff_locations [diff_locations(2:end)-1;end]])';
end
if (size(qrs_locations)>1)
    merge_ranges = find((qrs_locations(1,2:end)-qrs_locations(2,1:end-1)) < merge_window).*2;
    qrs_locations([merge_ranges merge_ranges+1]) = []; 
end
qrs_locations = reshape(qrs_locations,2,[])';
qrs_locations = round(mean(qrs_locations,2));

% Flag and drop low confidence detection ranges
valid_qrs = diff(qrs_locations)./movmean(diff(qrs_locations),5);
valid_qrs = [-1; find(valid_qrs<valid_threshold | valid_qrs>1/valid_threshold); length(valid_qrs)+1] + 1;
valid_ranges = [1 + valid_qrs(diff(valid_qrs)>=valid_window) ...
    valid_qrs(find(diff(valid_qrs)>=valid_window)+1) - 1];
valid_qrs_locations = cell2mat(arrayfun(@(x,y)x:y, ...
	valid_ranges(:,1),valid_ranges(:,2),'UniformOutput',0)');

valid_ranges = qrs_locations(valid_ranges);
qrs_locations = qrs_locations(valid_qrs_locations);

end
