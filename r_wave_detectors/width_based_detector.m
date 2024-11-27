function [ widths ] = width_based_detector(x,fs )

%% R-Wave filtering
[b,a] = butter(4, [5 45]/(fs/2), 'bandpass');
x=filtfilt(b,a,x);

clear b a

%% Initial peak detection

% We perform a static-threshold detection of R peaks during the first 10
% second on a half-wave rectified signal.

clamp=0;
x(find(x(1:10*fs)<clamp))=clamp;

clear clamp

initial_thr=7.5*mean(x((1:10*fs)))

for i=1:10*fs
    if x(i)>initial_thr
        thr_exceeded(i)=1;
    else
        thr_exceeded(i)=0;
    end
end

crossings=diff(thr_exceeded);
rising_crossings=find(crossings==1);
falling_crossings=find(crossings==-1);
all_peaks=ceil((rising_crossings+falling_crossings)/2);
widths=falling_crossings-rising_crossings;

%% Peak Detection
thr=initial_thr*ones(1,fs*10);
for i=(10*fs)+1:length(x)
    
    % Threshold Setting
     
    thr(i)=median(x(all_peaks(end-4:end)))/2; % The threshold is HALF of the median amplitude of the previous 5 peaks
       
    if x(i)>thr(i)
        thr_exceeded(i)=1;
    else
        thr_exceeded(i)=0;
    end
    
    % Checking if the threshold was crossed
    
    if thr_exceeded(i)~=thr_exceeded(i-1)
        
        if thr_exceeded(i)==1 % Rising crossing
            crossings(i)=1;
            rising_crossings(end+1)=i;
        end
        
        if thr_exceeded(i)==0 % Falling crossing - the peak is over.
            crossings(i)=-1;
            falling_crossings(end+1)=i;
            all_peaks(end+1)=ceil((rising_crossings(end)+falling_crossings(end))/2);
            widths(end+1)=falling_crossings(end) - rising_crossings(end);
                       
                
        end
        
    else % If the threshold was NOT crossed...
        
        crossings(i)=0;
    
    end
    
end

figure; plot(x); hold on; plot(thr); plot(all_peaks, x(all_peaks),'go'); 
figure; plot(widths)


end