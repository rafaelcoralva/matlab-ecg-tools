function [ r_peaks ] = r_peak_detector( ecg,original_sampling_frequency )
% R wave ECG peak detector. 

% First a band-pass filter is used to attenuate the entire signal outside of the R wave bandwidth. The NEO/TEO is
% then used to accentuate the R waves. A moving threshold that is esentially a smoothening of the rectification of 
% the TEO/NEO of the signal about the local median of the (normalised) TEO/NEO is then passed through the (normalised)
% TEO/NEO of the signal. The midpoint of the threshold crossings is then considered as the peak.

% Because this function uses the median to evaluate the threshold, it can be quite computationally costly on highly-sampled
% signals. Hence, if the sampling frequency is greater than 2kHz, the ECG
% is decimated by a factor of 5.

%% Parameters

window_width=2; % In seconds.
c=1.3; % Aggresivity factor.
refractory_period=0.2; % Corresponds to a maximum heart rate of 240 bpm - no R peaks can be detected in this window following the previous detection.


%% Decimation

if original_sampling_frequency>2000
    
    ecg=decimate(ecg,5);
    sampling_frequency=original_sampling_frequency/5;
    
end

if original_sampling_frequency<2001
    
    sampling_frequency=original_sampling_frequency;
    
end



%% Band-Pass Filtering

low_r_band=4;
high_r_band=30;

bpFilt = designfilt('bandpassiir','FilterOrder',10, ...
         'HalfPowerFrequency1',low_r_band,'HalfPowerFrequency2',high_r_band, ...
         'SampleRate',sampling_frequency);

filtered_ecg=filtfilt(bpFilt,ecg);

clear bpFilt


%% NEO of filtered signal

neo_ecg=neo(normalize(filtered_ecg)');
neo_ecg(end)=[]; % This is becasue the final element after NEOing is always zero, which can upset the normalization.
neo_ecg=normalize(neo_ecg);


%% Moving Threshold

hamm_window=0.2*sampling_frequency; % 60ms hamming window used to slide and smoothen through the rectified window signal.
hamm_filter=hamming(hamm_window);
hamm_filter=hamm_filter/sum(hamm_filter); % Normalizing the Hamming filter.

thr_exceeded=0*ones(1,length(neo_ecg)); % This vector holds the samples where the TEO of the filtered signal exceeds the moving threshold.

moving_median(1:window_width*sampling_frequency)=median(neo_ecg(1:window_width*sampling_frequency)); % The median for the first window_width seconds of the signal is static.

for i=window_width*sampling_frequency+1:1:length(neo_ecg) % The median for the rest of the signal is calculated in a lagging window of window_width seconds.
    
    moving_median(i)=median(neo_ecg(i-window_width*sampling_frequency:5:i));
    
end

thr=c*conv(hamm_filter,abs(neo_ecg-moving_median) + moving_median); % The threshold is essentially a smoothening of a the rectification of the TEO of the signal about its (moving) median.

% Adjusting the 'thr' vector following convolution to ensure it is
% synchronised to the signal, and of equal length.

thr(1:hamm_window/2)=[];
thr(end-hamm_window/2:end)=[];
thr(end-window_width/4:end)=thr(end-window_width);
thr(end:end+length(neo_ecg)-length(thr))=thr(end);

% Detecting the regions where the threshold is exceeded.

for i=1:length(neo_ecg)
    
    if neo_ecg(i)>thr(i)
        
        thr_exceeded(i)=1;
    end
    
end

% Interpolating between threshold crossing points

diff_thr_exceeded=diff(thr_exceeded); % This returns a vector of zeros except at the points where the signal CROSSED the threshold. When the threshold is first crossed (i.e. on the rising edge), diff_thr_exceeded is 1. Then it is zero until the smoothed 'x' crosses the threshold again on the falling edge, when it will become -1.

rising_edge_crossing=find(diff_thr_exceeded==1); % Find the rising edge threshold crossings
falling_edge_crossing=find(diff_thr_exceeded==-1); % Find the falling edge threshold crossings.

%% Checks

if rising_edge_crossing(1)>falling_edge_crossing(1) % If this is the case, the rising edge of the first peak has been missed because the signal starts too late. In this case, the first falling edge should be discarded.
    falling_edge_crossing(1)=[];
end

if rising_edge_crossing(end)>falling_edge_crossing(end) % If this is the case, the falling edge of the last peak has been missed becasuse the signal ended too soon. In this case, the last rising edge should be discarded.
    rising_edge_crossing(end)=[];
end

r_peaks=round((rising_edge_crossing+falling_edge_crossing)/2); % The wave peak position is the midpoint between falling and rising edge crossings.


% Refractory period check

rr_interval=diff(r_peaks);
r_peaks(find(rr_interval<refractory_period*sampling_frequency)+1)=[]; % Deleting the R peaks that occurred within the refractory period.

% Local maximum check - checks to see if the detected peak is the true peak
% (sometimes it occurs just before or just after the detected peak). This is
% a consequence of the interpolation (considering the peak as the midpoint
% between threshold crossings.

for i=1:length(r_peaks)
    
    if ecg(r_peaks(i)-1)>ecg(r_peaks(i))
        
        r_peaks(i)=r_peaks(i)-1;
    end
    
    if ecg(r_peaks(i)-2)>ecg(r_peaks(i))
        
        r_peaks(i)=r_peaks(i)-2;
    end
    
    if ecg(r_peaks(i)+1)>ecg(r_peaks(i))
        
        r_peaks(i)=r_peaks(i)+1;
    end
    
    if ecg(r_peaks(i)+2)>ecg(r_peaks(i))
        
        r_peaks(i)=r_peaks(i)+2;
    end
    
end

%% Plotting for troubleshooting

t=0:1/sampling_frequency:(length(neo_ecg)-1)/sampling_frequency;
r_peak_times=r_peaks/sampling_frequency;

% Filtered Signal

% figure
% plot(t,ecg(1:end-1))
% hold on
% plot(t,filtered_ecg(1:end-1))
% xlabel('Time (s)')
% ylabel('mV')
% legend('Original Signal','Bandpass-Filtered Signal')
% title('Original vs. Filtered Signal')

% TEO of the Signal

figure
plot(t,neo_ecg)
hold on
plot(t,thr)
plot(r_peak_times,neo_ecg(r_peaks),'o')
xlabel('Time (s)')
ylabel('Arb Units')
legend('Normalised NEO of BP ECG','Moving-median Threshold')
title('Normalised NEO/TEO of the band-pass filtered ECG')

% figure
% plot(t,ecg(1:end-1))
% hold on
% plot(r_peak_times,ecg(r_peaks),'o')
% xlabel('Time (s)')
% ylabel('Amplitude (mV)')
% title('R peak-labelled ECG')


if original_sampling_frequency>2000
    
    r_peaks=r_peaks*5; % This is necessary to adjust for the decimation.
    
end

end

