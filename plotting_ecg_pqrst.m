function plotting_ecg_pqrst(coherent_ecg_pqrst, preprocessed_ecg, fs)
% Inside of Wrapped Simulink Function to plot the PQRST markers of the
% preprocessed ECG.

t = 0:1/fs:((length(preprocessed_ecg)-1))/fs;

coherent_ecg_P = coherent_ecg_pqrst(isfinite(coherent_ecg_pqrst(:,1)),1);
coherent_ecg_Q = coherent_ecg_pqrst(isfinite(coherent_ecg_pqrst(:,2)),2);
coherent_ecg_R = coherent_ecg_pqrst(isfinite(coherent_ecg_pqrst(:,3)),3);
coherent_ecg_S = coherent_ecg_pqrst(isfinite(coherent_ecg_pqrst(:,4)),4);
coherent_ecg_T = coherent_ecg_pqrst(isfinite(coherent_ecg_pqrst(:,5)),5);

fig_ecg_labelled = figure;
plot(t,preprocessed_ecg,'b'); title('Preprocessed ECG'); xlabel('Time (sec)'); ylabel('Amplitude (mV)'); box off; scrollplot; hold on;
plot(coherent_ecg_P/fs,preprocessed_ecg(coherent_ecg_P),'ro','LineWidth',1.5);
plot(coherent_ecg_Q/fs,preprocessed_ecg(coherent_ecg_Q),'o','LineWidth', 1.5,'color',[1 0.75 0]);
plot(coherent_ecg_R/fs,preprocessed_ecg(coherent_ecg_R),'go','LineWidth', 1.5);
plot(coherent_ecg_S/fs,preprocessed_ecg(coherent_ecg_S),'bo','LineWidth', 1.5);
plot(coherent_ecg_T/fs,preprocessed_ecg(coherent_ecg_T),'ko','LineWidth', 1.5);

% Plotting the legend (depends on the which peaks were detected).
if ~isempty(coherent_ecg_P) & ~isempty(coherent_ecg_Q) & ~isempty(coherent_ecg_R) & ~isempty(coherent_ecg_S) & ~isempty(coherent_ecg_T)
    legend('ecg','p','q','r','s','t')
elseif isempty(coherent_ecg_P) & ~isempty(coherent_ecg_Q) & ~isempty(coherent_ecg_R) & ~isempty(coherent_ecg_S) & ~isempty(coherent_ecg_T) 
    legend('ecg','q','r','s','t')  
elseif isempty(coherent_ecg_P) & isempty(coherent_ecg_Q) & ~isempty(coherent_ecg_R) & isempty(coherent_ecg_S) & isempty(coherent_ecg_T) 
    legend('ecg','r')
end

end

