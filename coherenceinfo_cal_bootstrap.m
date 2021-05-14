function Coherence = coherenceinfo_cal_bootstrap(CoherencyT_unshifted,Freqs, CSR, CSR_low, Nyquist,Fs,nFFT,TR,FeatureName)
%function Coherence = coherenceinfo_cal(CoherencyT_unshifted,Freqs, CSR, CSR_up, CSR_low, Nyquist,Fs,nFFT,TR,PlotCoherenceFig, FeatureName)

% PlotCoherenceFig=0;

if nargin<10
    FeatureName = 'Feature';
end

CoherencyT = fftshift(CoherencyT_unshifted);

% Calculate the information on coherence

% normalize coherencies
Freqs4Info = Freqs;
Coherence_value = CSR.^2;
Coherence4Info = Coherence_value;
% Coherence_up = CSR_up.^2;
% CoherenceUp4Info = Coherence_up;

closgn = sign(real(CSR_low));
Coherence_low = (CSR_low.^2) .* closgn; %Coherence_low{cc} can be negative, multiply by sign after squaring
CoherenceLow4Info = Coherence_low;

%% restrict frequencies analyzed to the requested cutoff and minimum frequency given the window size
freqCutoff = Nyquist; % keep all frequencies up to Nyquist as of now
if freqCutoff ~= -1
    eindx = find(Freqs4Info < freqCutoff,1,'last');
    indx = 1:eindx;
    
    Freqs4Info = Freqs4Info(indx);
    Coherence4Info = Coherence4Info(indx);
%     CoherenceUp4Info = CoherenceUp4Info(indx);
    CoherenceLow4Info = CoherenceLow4Info(indx);
end

% The DC component is now eliminated by zscoring we can discard that
% section
% minFreq = Fs/nFFT;
% if minFreq > 0
%     findx = find(Freqs4Info >= minFreq); %#ok<MXFND>
%     sindx = min(findx);
%     Freqs4Info = Freqs4Info(sindx:end);
%     Coherence4Info = Coherence4Info(sindx:end);
% %     CoherenceUp4Info = CoherenceUp4Info(sindx:end);
%     CoherenceLow4Info = CoherenceLow4Info(sindx:end);
% end

%% if the clower goes below zero set all values to zero and keep track fo that first value of frequency for which the coherence is non significant
cutoffIndex = find(CoherenceLow4Info < 0, 1, 'first');
if (isempty(cutoffIndex))
    cutoffIndex = length(CoherenceLow4Info);
elseif cutoffIndex==1 % The first frequency is not significant, try to find the next non significant that follows the significant ones
    cutoffIndex = find(find(CoherenceLow4Info < 0)> find(CoherenceLow4Info > 0,1,'first'),1,'first');
end
freqCutoff = Freqs4Info(cutoffIndex);
FirstNonSigCoherenceFreq = freqCutoff;

%% Calculate the coherence weighted significant frequencies
Weight = Coherence4Info(CoherenceLow4Info>0);
Weight = Weight./sum(Weight);
Freqs4weight = Freqs4Info(CoherenceLow4Info>0);
WeightedSigCoherenceFreq = sum(Freqs4weight.*Weight);

%% Calculate the cumulative sum of the significant information pdf
CumSumWeight = cumsum(Weight); 
CumSumSigCoherence50Hz = CumSumWeight(find(Freqs4weight<=50,1,'Last'));

%% compute information by integrating log of 1 - coherence
df = Freqs4Info(2) - Freqs4Info(1);
Info = -df*sum(log2(1 - Coherence4Info(CoherenceLow4Info>0)));
% Info_up = -df*sum(log2(1 - CoherenceUp4Info(1:cutoffIndex)));
% Info_low = -df*sum(log2(1 - CoherenceLow4Info(1:cutoffIndex)));

%% Lowpass the coherency according to threshold FirstNonSigCoherenceFreq and find the width of max peak
if FirstNonSigCoherenceFreq
    [z,p,k] = butter(6,FirstNonSigCoherenceFreq/(Fs/2),'low');
    sos_low = zp2sos(z,p,k);
    CoherencyT_filt=filtfilt(sos_low,1,CoherencyT);
    CoherencyT_xTimeDelay = -(((nFFT/2)/Fs)*10^3):TR:(((nFFT/2-1))/Fs)*10^3; % Corresponding values in ms of the Delay for each value of CoherencyT
    [P,Locs] = findpeaks(CoherencyT_filt);

    if ~isempty(Locs)
    %     [P,IndM] = max(P);
        [~,IndM] = max(P);
        Locs = Locs(IndM);
        CoherencyT_DelayAtzero = CoherencyT_xTimeDelay(Locs);
        CrossZero1 = CoherencyT_xTimeDelay(find(CoherencyT_filt(1:Locs)<=mean(CoherencyT_filt), 1, 'last')+1);
        CrossZero2 = CoherencyT_xTimeDelay(Locs + find(CoherencyT_filt(Locs+1:end)<=mean(CoherencyT_filt), 1, 'first')-1);
        if isempty(CrossZero2)
            CrossZero2 = CoherencyT_xTimeDelay(Locs + find(CoherencyT_filt(Locs+1:end)==min(CoherencyT_filt(Locs+1:end)), 1, 'first')-1);
        end
        CoherencyT_WidthAtMaxPeak = CrossZero2 - CrossZero1;
    else
        CoherencyT_DelayAtzero = nan;
        CoherencyT_WidthAtMaxPeak = nan;
        warning('There is no peak in Coherency T \n')
    end
else
    CoherencyT_DelayAtzero = nan;
        CoherencyT_WidthAtMaxPeak = nan;
        warning('There is no peak in Coherency T \n')
end

%% Plot the value of coherency as a function of delay
% if PlotCoherenceFig
%     figure(3)
%     clf
%     subplot(1,3,1)
%     plot(CoherencyT_xTimeDelay, CoherencyT_filt, 'LineWidth',2);
%     if ~isempty(P)
%         hold on
%         plot(CoherencyT_DelayAtzero,P,'ro', 'MarkerSize',16, 'MarkerFaceColor','r')
%         text(max(CoherencyT_xTimeDelay)/5, 4.5/5*P, sprintf('Delay at zero = %.2f ms', CoherencyT_DelayAtzero));
%         text(max(CoherencyT_xTimeDelay)/5, 4/5*P, sprintf('Width at y=mean(coherencyT) = %.2f ms', CoherencyT_WidthAtMaxPeak));
%         hold off
%     end
%     xlabel('Time Delay')
%     ylabel(sprintf('Coherency with %s', FeatureName))
%     %     subplot(1,2,2)
%     %     plot(Freqs, Coherence_value, '-k','LineWidth',2)
%     %     hold on
%     %     plot(Freqs, Coherence_up, '--r','LineWidth',2)
%     %     hold on
%     %     plot(Freqs, Coherence_low, '--b','LineWidth',2)
%     %     hold off
%     subplot(1,3,2)
%     shadedErrorBar(Freqs, Coherence_value,[(Coherence_up-Coherence_value)'; (-Coherence_low+Coherence_value)'], {'LineWidth',2,'Color','k'});
%     hold on
%     plot([Freqs(1) Freqs(end)], [0 0], 'r--', 'LineWidth',2);
%     MaxY = 0.5;
%     ylim([-0.1 MaxY]);
%     text(Freqs(end)/3,0.8*MaxY,sprintf('Info = %.2f   InfoUp = %.2f    InfoLow = %.2f', Info, Info_up, Info_low));
%     hold on
%     text(Freqs(end)/3,0.85*MaxY,sprintf('FreqCutOff = %.2f Hz -> TimeResolution = %.2f ms', FirstNonSigCoherenceFreq,1/FirstNonSigCoherenceFreq*10^3));
%     hold off
%     xlabel('Frequencies (Hz)');
%     ylabel(sprintf('Coherence with %s', FeatureName));
%     
%     
%     if ~isempty(FirstNonSigCoherenceFreq)
%         if ~isempty(CoherencePeaksF)
%             maxFreqInd = max(find(Freqs == max(CoherencePeaksF)), find(FirstNonSigCoherenceFreq==Freqs)) +2;
%         else
%             maxFreqInd = find(FirstNonSigCoherenceFreq==Freqs) +2;
%         end
%         subplot(1,3,3)
%         shadedErrorBar(Freqs(1:maxFreqInd), Coherence_value(1:maxFreqInd),[(Coherence_up(1:maxFreqInd)-Coherence_value(1:maxFreqInd))'; (-Coherence_low(1:maxFreqInd)+Coherence_value(1:maxFreqInd))'], {'LineWidth',2,'Color','k'});
%         hold on
%         plot([Freqs(1) Freqs(maxFreqInd)], [0 0], 'r--', 'LineWidth',2);
%         hold on
%         plot(CoherencePeaksF, CoherencePeaks, 'go', 'MarkerSize',6,'MarkerFaceColor','g');
%         hold off
%         xlabel('Frequencies (Hz)');
%         ylabel(sprintf('Coherence with %s', FeatureName));
%     end
% end

Coherence.Info = Info;
Coherence.CoherencyT_DelayAtzero = CoherencyT_DelayAtzero;
Coherence.CoherencyT_WidthAtMaxPeak = CoherencyT_WidthAtMaxPeak;
Coherence.WeightedSigCoherenceFreq = WeightedSigCoherenceFreq;
Coherence.CumSumSigCoherence50Hz = CumSumSigCoherence50Hz;
Coherence.FeatureName = FeatureName;

end