function [t_peaks ] = negative_t_peak_detector_that_depends_on_r_peaks( ecg,r_peaks,sampling_frequency)
% Funciton to detect the T wave of the ECG based on the already detected R wave peak positions.

% This function looks for the T wave in a window after the R peak. This window spans 10-50% of the previous RR interval, following each detected R peak,
% scaled by the duration of the previous RR interval. These windows are normalized and static thresholds are set at 0.95. The peaks are then calculated by
% linearly interpolating between rising and falling threshold crossings.

% Negative T wave peaks are accounted for by firstly rectifying the ECG.

%% Parameters

lower_window_limit=0.1;
higher_window_limit=0.5;

thr=0.95;
rr_interval=diff(r_peaks);


%% Band-Pass Filtering

low_r_band=0.5;
high_r_band=20;

bpFilt = designfilt('bandpassiir','FilterOrder',10, ...
         'HalfPowerFrequency1',low_r_band,'HalfPowerFrequency2',high_r_band, ...
         'SampleRate',sampling_frequency);

filtered_ecg=filtfilt(bpFilt,ecg);

clear bpFilt

%% First Beat Detection 

% For the fist T wave peak we consider the following RR interval (because we don't know it's preceding RR interval).

t_window=filtered_ecg((r_peaks(1)+(rr_interval(1)*lower_window_limit)):(r_peaks(1)+(rr_interval(1)*higher_window_limit)));
t_window=abs(t_window-mean(t_window)); % Rectification is performed on a zero-mean window. The ECG is generally zero-mean by definition but this ensures any wandering baseline doesn't upset any peak detection following signal rectification.


norm_t_window=normalize(t_window);

thr_exceeded=0*ones(1,length(t_window));

for j=1:length(norm_t_window)
    
    if norm_t_window(j)>thr
        thr_exceeded(j)=1;
    end
    
end

diff_thr_exceeded=diff(thr_exceeded); % This returns a vector of zeros except at the points where the signal CROSSED the threshold. When the threshold is first crossed (i.e. on the rising edge), diff_thr_exceeded is 1. Then it is zero until the signal crosses the threshold again on the falling edge, when it will become -1.

rising_edge_crossing=find(diff_thr_exceeded==1); % Find the rising edge threshold crossings
falling_edge_crossing=find(diff_thr_exceeded==-1); % Find the falling edge threshold crossings.


% Check

peak=peak_checker(thr,rising_edge_crossing,falling_edge_crossing,norm_t_window);

t_peaks=r_peaks(1)+(lower_window_limit*rr_interval(1)) + peak; % The (1) index for the rising and falling edge crossings accounts for cases where there are several "peaks" or when the repolarization baseline is quite elevated, and acts as the "peak".

%% Detection

% For the following T wave peaks we consider the preceding RR interval when
% opening a window where we look for the T peak.

for i=2:length(r_peaks) % We cannot open the T wave window of the first beat because we dont know its preceding RR interval.
    
    t_window=filtered_ecg((r_peaks(i)+rr_interval(i-1)*lower_window_limit):(r_peaks(i)+rr_interval(i-1)*higher_window_limit));
    t_window=abs(t_window-mean(t_window)); % Rectification is performed on a zero-mean window. The ECG is generally zero-mean by definition but this ensures any wandering baseline doesn't upset any peak detection following signal rectification.
    
    norm_t_window=normalize(t_window);
    
    thr_exceeded=0*ones(1,length(t_window));
    
    for j=1:length(t_window)
        
        if norm_t_window(j)>thr
            thr_exceeded(j)=1;
        end
        
    end
    
    diff_thr_exceeded=diff(thr_exceeded); % This returns a vector of zeros except at the points where the signal CROSSED the threshold. When the threshold is first crossed (i.e. on the rising edge), diff_thr_exceeded is 1. Then it is zero until the signal crosses the threshold again on the falling edge, when it will become -1.
    
    rising_edge_crossing=find(diff_thr_exceeded==1); % Find the rising edge threshold crossings.
    falling_edge_crossing=find(diff_thr_exceeded==-1); % Find the falling edge threshold crossings.
    
    % Checking the peak detection
    
    peak=peak_checker(thr,rising_edge_crossing,falling_edge_crossing,norm_t_window);
    
    t_peaks(i)=r_peaks(i)+(lower_window_limit*rr_interval(i-1)) + peak; 
    
end

t_peaks=ceil(t_peaks);

%% Plotting for troubleshooting

% Plotting the effect of Filtering

t=0:1/sampling_frequency:(length(ecg)-1)/sampling_frequency;

figure;
plot(t,ecg)
hold on
plot(t,filtered_ecg)
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Raw and Filtered ECG Signal for T Wave detection')

% sampling_frequency=10000; % CHECK THIS IS TRUE - THE FUNCTION DOESN'T TAKE THE SAMPLING FREQUENCY AS INPUT.

% Plotting the detected T-waves on the filtered signal

t_peak_times=t_peaks/sampling_frequency;

figure
plot(t,filtered_ecg)
hold on
plot(t_peak_times,filtered_ecg(t_peaks),'o')
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('ECG with labelled T wave peaks')

end

