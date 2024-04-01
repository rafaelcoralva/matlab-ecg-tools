function [p_peaks ] = p_peak_detector_that_depends_on_r_peaks( ecg,r_peaks,sampling_frequency)
% Funciton to detect the P wave of the ECG based on the already detected R wave peak positions.

% This function looks for the P wave in a window just before the R peak. This window spans 10-35% of the previous RR interval, preceding each detected R peak,
% scaled by the duration of the previous RR interval. First, the median of these windows is calculated to see if the P wave is positive or negative. If it is
% negative, the window is flipped (x-1). These windows are then normalized and static thresholds are set at 0.98. The peaks are then calculated by
% linearly interpolating between rising and falling threshold crossings.



%% Parameters

higher_window_limit=0.35;
lower_window_limit=0.1;

thr=0.98;
rr_interval=diff(r_peaks);


%% Band-Pass Filtering

low_r_band=1;
high_r_band=30;

bpFilt = designfilt('bandpassiir','FilterOrder',6, ...
    'HalfPowerFrequency1',low_r_band,'HalfPowerFrequency2',high_r_band, ...
    'SampleRate',sampling_frequency);

filtered_ecg=filtfilt(bpFilt,ecg);

clear bpFilt

%% First Beat Detection

% For the fist P wave peak we consider the following RR interval (because we don't know it's preceding RR interval).

if (r_peaks(1)-(rr_interval(1)*higher_window_limit))<0
    
    r_peaks(1)=[];
    
end

p_window=filtered_ecg((r_peaks(1)-(rr_interval(1)*higher_window_limit)):(r_peaks(1)-(rr_interval(1)*lower_window_limit)));

if median(p_window)<0
    p_window=p_window*-1;
end

norm_p_window=normalize(p_window);

thr_exceeded=0*ones(1,length(p_window));

for j=1:length(norm_p_window)
    
    if norm_p_window(j)>thr
        thr_exceeded(j)=1;
    end
    
end

diff_thr_exceeded=diff(thr_exceeded); % This returns a vector of zeros except at the points where the signal CROSSED the threshold. When the threshold is first crossed (i.e. on the rising edge), diff_thr_exceeded is 1. Then it is zero until the signal crosses the threshold again on the falling edge, when it will become -1.

rising_edge_crossing=find(diff_thr_exceeded==1); % Find the rising edge threshold crossings
falling_edge_crossing=find(diff_thr_exceeded==-1); % Find the falling edge threshold crossings.


% Check

peak=peak_checker(thr,rising_edge_crossing,falling_edge_crossing,norm_p_window);

p_peaks=ceil(r_peaks(1)-(higher_window_limit*rr_interval(1)) + peak); % The (1) index for the rising and falling edge crossings accounts for cases where there are several "peaks" or when the repolarization baseline is quite elevated, and acts as the "peak".

%% Detection

% For the following T wave peaks we consider the preceding RR interval when
% opening a window where we look for the T peak.

for i=2:length(r_peaks) % We cannot open the T wave window of the first beat because we dont know its preceding RR interval.
    
    p_window=filtered_ecg((r_peaks(i)-rr_interval(i-1)*higher_window_limit):(r_peaks(i)-rr_interval(i-1)*lower_window_limit));
    
    if median(p_window)<0
        p_window=p_window*-1;
    end
    
    norm_p_window=normalize(p_window);
    
    thr_exceeded=0*ones(1,length(p_window));
    
    for j=1:length(p_window)
        
        if norm_p_window(j)>thr
            thr_exceeded(j)=1;
        end
        
    end
    
    diff_thr_exceeded=diff(thr_exceeded); % This returns a vector of zeros except at the points where the signal CROSSED the threshold. When the threshold is first crossed (i.e. on the rising edge), diff_thr_exceeded is 1. Then it is zero until the signal crosses the threshold again on the falling edge, when it will become -1.
    
    rising_edge_crossing=find(diff_thr_exceeded==1); % Find the rising edge threshold crossings.
    falling_edge_crossing=find(diff_thr_exceeded==-1); % Find the falling edge threshold crossings.
    
    % Checking the peak detection
    
    peak=peak_checker(thr,rising_edge_crossing,falling_edge_crossing,norm_p_window);
    
    p_peaks(i)=ceil(r_peaks(i)-(higher_window_limit*rr_interval(i-1)) + peak);
    
    
    % Local Maxima Check
    
    if ecg(p_peaks(i)-1)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)-1;
    end
    
    if ecg(p_peaks(i)-2)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)-2;
    end
    
    if ecg(p_peaks(i)-3)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)-3;
    end
    
    if ecg(p_peaks(i)-4)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)-4;
    end
    
    if ecg(p_peaks(i)-5)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)-5;
    end
    
    
    if ecg(p_peaks(i)+1)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)+1;
    end
    
    if ecg(p_peaks(i)+2)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)+2;
    end
    
    if ecg(p_peaks(i)+3)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)+3;
    end
    
    if ecg(p_peaks(i)+4)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)+4;
    end
    
    if ecg(p_peaks(i)+5)>ecg(p_peaks(i))
        p_peaks(i)=p_peaks(i)+5;
    end
    
    
end


%% Plotting for troubleshooting

% Plotting the effect of Filtering
%
% t=0:1/sampling_frequency:(length(ecg)-1)/sampling_frequency;
%
% figure;
% plot(t,ecg)
% hold on
% plot(t,filtered_ecg)
% xlabel('Time (s)')
% ylabel('Amplitude (mV)')
% title('Raw and Filtered ECG Signal for P Wave detection')

% Plotting the detected T-waves on the filtered signal
%
% t_peak_times=t_peaks/sampling_frequency;
%
% figure
% plot(t,filtered_ecg)
% hold on
% plot(p_peak_times,filtered_ecg(p_peaks),'o')
% xlabel('Time (s)')
% ylabel('Amplitude (mV)')
% title('ECG with labelled P wave peaks')

end

