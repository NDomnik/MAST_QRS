% MAST_QRS: Automated detection of ECG R-waves, optimized for murine ECG 
% obtained via implantable radiotelemetry
% Copyright (C) 2020 Sami Torbey

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

function [ qrs ] = mast_qrs( ecg, freq )

if (nargin<2 || isempty(freq))
    freq = 1000;
else
    freq = freq(1);
end

if (size(ecg,2) > size(ecg,1))
    ecg = ecg';
end

len_ecg = size(ecg,1)-1;
num_leads = size(ecg,2);

narrow_window = round(0.2*freq);
wide_window = round(3*freq);
merge_window = round(0.05*freq);
discard_window = round(0.05*freq);
decline_rate = 20/freq;

[S_numer,S_denom] = butter(4,[7 17]*2/freq);
ecg = filtfilt(S_numer,S_denom,ecg);
ecg = abs(diff(ecg));

narrow_threshold = fast_smooth(ecg,narrow_window);
wide_threshold = (fast_smooth(ecg,wide_window) + ...    
    fast_smooth(ecg,round(sqrt(narrow_window*wide_window))) + ...
    narrow_threshold) / 3;

rr_coeff = zeros(size(ecg));
for i=1:num_leads
    qrs = find_ranges(narrow_threshold(:,i)>wide_threshold(:,i),merge_window);
    rr = diff([0;qrs(:,1)]);
    rr = abs(rr./fast_smooth(rr,5)-1);
    rr = slidefun(@max,5,10.^rr-1);
    rr_coeff(:,i) = interp1([0;mean(qrs,2);len_ecg+1],rr([1 1:end end]),1:len_ecg,'nearest')';    
    
    for j=1:size(qrs,1)
        wide_threshold(qrs(j,1):qrs(j,2),i) = max(wide_threshold(qrs(j,1):qrs(j,2),i));
    end
    for j=2:len_ecg
        if narrow_threshold(j,i)<=wide_threshold(j-1,i)
            wide_threshold(j,i)=decline_rate*wide_threshold(j,i)+(1-decline_rate)*wide_threshold(j-1,i);
        end
    end
end

narrow_threshold = (narrow_threshold+rr_coeff.*wide_threshold) ./ (rr_coeff+1);
narrow_threshold = mean((narrow_threshold-wide_threshold)./(narrow_threshold+wide_threshold),2);
narrow_threshold = fast_smooth(narrow_threshold,narrow_window);
narrow_threshold = narrow_threshold .* min(rr_coeff+1,[],2);

qrs = find_ranges(narrow_threshold>0 & moving_std(narrow_threshold,wide_window)>0.003,merge_window,discard_window);
qrs = round(mean(qrs,2));

end
