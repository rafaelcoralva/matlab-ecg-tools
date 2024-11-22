function plotting_ecg_pqrst(idxs_coherent_pqrst, preprocessed_ecg, fs)
% Wrapped Simulink Function to plot the PQRST markers of the preprocessed ECG.
% INPUTS:
% idxs_coherent_pqrst (matrix: [Nx5]) - Matrix where each column represents P,Q,R,S,T index positions in the ECG signal.
% preprocessed_ecg (vector) - Preprocessed ECG signal to be plotted
% fs (positive integer scalar) - ECG sampling frequency in Hz.

% Creating time vector.
t = 0:1/fs:((length(preprocessed_ecg)-1))/fs;

% Extracting indexes of P,Q,R,S,T waves/
idxs_coherent_P = idxs_coherent_pqrst(isfinite(idxs_coherent_pqrst(:,1)),1);
idxs_coherent_Q = idxs_coherent_pqrst(isfinite(idxs_coherent_pqrst(:,2)),2);
idxs_coherent_R = idxs_coherent_pqrst(isfinite(idxs_coherent_pqrst(:,3)),3);
idxs_coherent_S = idxs_coherent_pqrst(isfinite(idxs_coherent_pqrst(:,4)),4);
idxs_coherent_T = idxs_coherent_pqrst(isfinite(idxs_coherent_pqrst(:,5)),5);

% Plotting labelled ECG.
fig_ecg_labelled = figure;
plot(t,preprocessed_ecg,'b'); title('Preprocessed ECG'); xlabel('Time (sec)'); ylabel('Amplitude (mV)'); box off; scrollplot; hold on;
plot(idxs_coherent_P/fs,preprocessed_ecg(idxs_coherent_P),'ro','LineWidth', 1.5);
plot(idxs_coherent_Q/fs,preprocessed_ecg(idxs_coherent_Q),'o','LineWidth', 1.5,'color',[1 0.75 0]);
plot(idxs_coherent_R/fs,preprocessed_ecg(idxs_coherent_R),'go','LineWidth', 1.5);
plot(idxs_coherent_S/fs,preprocessed_ecg(idxs_coherent_S),'bo','LineWidth', 1.5);
plot(idxs_coherent_T/fs,preprocessed_ecg(idxs_coherent_T),'ko','LineWidth', 1.5);

% Plotting the legend (depends on the which peaks were detected).
if ~isempty(idxs_coherent_P) & ~isempty(idxs_coherent_Q) & ~isempty(idxs_coherent_R) & ~isempty(idxs_coherent_S) & ~isempty(idxs_coherent_T)
    legend('ecg','p','q','r','s','t')
elseif isempty(idxs_coherent_P) & ~isempty(idxs_coherent_Q) & ~isempty(idxs_coherent_R) & ~isempty(idxs_coherent_S) & ~isempty(idxs_coherent_T) 
    legend('ecg','q','r','s','t')  
elseif isempty(idxs_coherent_P) & isempty(idxs_coherent_Q) & ~isempty(idxs_coherent_R) & isempty(idxs_coherent_S) & isempty(idxs_coherent_T) 
    legend('ecg','r')
end

end

