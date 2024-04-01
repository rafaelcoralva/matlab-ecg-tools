function [ peak_to_peak_amplitude ] = ecg_peak_to_peak_amplitude_calculator( sonr, peaks, troughs )
% Simple algorithm to calculate the peak-to-peak amplitude of a SonR
% component based on the input peak and trough times.

%   Peak-to-peak amplitude is simply peaks-troughs. This function assumes
%   that the index position of peaks and troughs are synchronized.

peak_to_peak_amplitude=sonr(peaks)-sonr(troughs);


end

