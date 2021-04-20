function Coherence = coherenceinfo_cal(CoherencyT_unshifted,Freqs, CSR, CSR_up, CSR_low, Nyquist,Fs,nFFT,TR,PlotCoherenceFig, FeatureName)
if nargin<10
    PlotCoherenceFig=0;
end
if nargin<11
    FeatureName = 'Feature';
end
MinCoherence4peaks = 0.01;

CoherencyT = fftshift(CoherencyT_unshifted);

% Calculate the information on coherence

% normalize coherencies
Freqs4Info = Freqs;
Coherence_value = CSR.^2;
Coherence4Info = Coherence_value;
Coherence_up = CSR_up.^2;
CoherenceUp4Info = Coherence_up;

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
    CoherenceUp4Info = CoherenceUp4Info(indx);
    CoherenceLow4Info = CoherenceLow4Info(indx);
end

minFreq = Fs/nFFT;
if minFreq > 0
    findx = find(Freqs4Info >= minFreq); %#ok<MXFND>
    sindx = min(findx);
    Freqs4Info = Freqs4Info(sindx:end);
    Coherence4Info = Coherence4Info(sindx:end);
    CoherenceUp4Info = CoherenceUp4Info(sindx:end);
    CoherenceLow4Info = CoherenceLow4Info(sindx:end);
end

%% if the clower goes below zero set all values to zero and keep track fo that first value of frequency for which the coherence is non significant
cutoffIndex = find(CoherenceLow4Info < 0, 1, 'first');
if (isempty(cutoffIndex))
    cutoffIndex = length(CoherenceLow4Info);
end
freqCutoff = Freqs4Info(cutoffIndex);
FirstNonSigCoherenceFreq = freqCutoff;

%% compute information by integrating log of 1 - coherence
df = Freqs4Info(2) - Freqs4Info(1);
Info = -df*sum(log2(1 - Coherence4Info(1:cutoffIndex)));
Info_up = -df*sum(log2(1 - CoherenceUp4Info(1:cutoffIndex)));
Info_low = -df*sum(log2(1 - CoherenceLow4Info(1:cutoffIndex)));

%% Lowpass the coherency according to threshold FirstNonSigCoherenceFreq and find the width of max peak
[z,p,k] = butter(6,FirstNonSigCoherenceFreq/(Fs/2),'low');
sos_low = zp2sos(z,p,k);
CoherencyT_filt=filtfilt(sos_low,1,CoherencyT);
CoherencyT_xTimeDelay = -(((nFFT/2)/Fs)*10^3):TR:(((nFFT/2-1))/Fs)*10^3; % Corresponding values in ms of the Delay for each value of CoherencyT
[P,Locs] = findpeaks(CoherencyT_filt);

if ~isempty(Locs)
    [P,IndM] = max(P);
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

%% Find if there are other peaks in coherence than the initial one
MaxCoherence1 = max(Coherence_value);
MaxCoherence2 = Freqs(Coherence_value==MaxCoherence1);
[Coherence2ndPeaks, LocsC] = findpeaks(Coherence_value, 'MinPeakHeight', MinCoherence4peaks, 'MinPeakProminence',MinCoherence4peaks);
LocsC(Coherence2ndPeaks==MaxCoherence1)=[];% eliminate the peak that corresponds to max value
Coherence2ndPeaks(Coherence2ndPeaks==MaxCoherence1)=[]; % eliminate the peak that corresponds to max value
Coherence2ndPeaks(Coherence_low(LocsC)<0)=[]; % eliminate peaks that are non-significant
LocsC(Coherence_low(LocsC)<0)=[];% eliminate peaks that are non-significant
CoherencePeaks = Coherence2ndPeaks';
CoherencePeaksF = Freqs(LocsC)';
if ~isempty(CoherencePeaksF)
    SecondCoherenceFreqCutOff = max(CoherencePeaksF);
else
    SecondCoherenceFreqCutOff = nan;
end


%% That was my own calculation without multitaper and jackknife
%     % Calculate cross-correlation between sound feature and neural data along with autocorrelation of neural data and sound feature for each stim
%     AcorrAmp = nan(NStims,length(Lags));
%     AcorrY = nan(NStims,length(Lags));
%     XcorrAmpY = nan(NStims,length(Lags));
%     W =hann(length(Lags))';
%     for ss=1:NStims
%         AcorrAmp(ss,:) = W.*xcorr(XAmpPerStim{ss},XAmpPerStim{ss},Delay,'unbiased');
%         AcorrY(ss,:) = W.*xcorr(YPerStim{ss},YPerStim{ss},Delay,'unbiased');
%         XcorrAmpY(ss,:) = W.*xcorr(XAmpPerStim{ss},YPerStim{ss},Delay,'unbiased');
%     end
%
%     % Getting the fft of the mean over stims for each crosscorrelation
%     FFT_AcorrAmp = fft(mean(AcorrAmp));
%     FFT_AcorrY = fft(mean(AcorrY));
%     FFT_XcorrAmpY = fft(mean(XcorrAmpY));
%
%     % calculate the coherency as a function of frequencies. In this domain,
%     % the values of power at each frequency band are independant
%     CoherencyF = (FFT_XcorrAmpY)./(abs(FFT_AcorrAmp) .* abs(FFT_AcorrY)).^.5; % here we take abs(FFT) for autocorrelation to alleviate the effect of the phase (pick of the autocorrelation should be at zero phase and it's not when calculating the fft)
%     Coherence_value = (abs(CoherencyF{cc})).^2;
%
%     % revert to time domain to find the coherency
%     CoherencyT = ifft(CoherencyF);

% Plot the value of coherency as a function of delay
if PlotCoherenceFig
    figure(3)
    clf
    subplot(1,3,1)
    plot(CoherencyT_xTimeDelay, CoherencyT_filt, 'LineWidth',2);
    if ~isempty(P)
        hold on
        plot(CoherencyT_DelayAtzero,P,'ro', 'MarkerSize',16, 'MarkerFaceColor','r')
        text(max(CoherencyT_xTimeDelay)/5, 4.5/5*P, sprintf('Delay at zero = %.2f ms', CoherencyT_DelayAtzero));
        text(max(CoherencyT_xTimeDelay)/5, 4/5*P, sprintf('Width at y=mean(coherencyT) = %.2f ms', CoherencyT_WidthAtMaxPeak));
        hold off
    end
    xlabel('Time Delay')
    ylabel(sprintf('Coherency with %s', FeatureName))
    %     subplot(1,2,2)
    %     plot(Freqs, Coherence_value, '-k','LineWidth',2)
    %     hold on
    %     plot(Freqs, Coherence_up, '--r','LineWidth',2)
    %     hold on
    %     plot(Freqs, Coherence_low, '--b','LineWidth',2)
    %     hold off
    subplot(1,3,2)
    shadedErrorBar(Freqs, Coherence_value,[(Coherence_up-Coherence_value)'; (-Coherence_low+Coherence_value)'], {'LineWidth',2,'Color','k'});
    hold on
    plot([Freqs(1) Freqs(end)], [0 0], 'r--', 'LineWidth',2);
    hold on
    plot(CoherencePeaksF, CoherencePeaks, 'go', 'MarkerSize',6,'MarkerFaceColor','g');
    hold on
    plot(MaxCoherence2, MaxCoherence1, 'ro', 'MarkerSize',6,'MarkerFaceColor','r');
    MaxY = 0.5;
    ylim([-0.1 MaxY]);
    text(Freqs(end)/3,0.8*MaxY,sprintf('Info = %.2f   InfoUp = %.2f    InfoLow = %.2f', Info, Info_up, Info_low));
    hold on
    text(Freqs(end)/3,0.9*MaxY,sprintf('MaxCoherence = %.2f', MaxCoherence1));
    hold on
    text(Freqs(end)/3,0.85*MaxY,sprintf('FreqCutOff = %.2f Hz -> TimeResolution = %.2f ms', FirstNonSigCoherenceFreq,1/FirstNonSigCoherenceFreq*10^3));
    hold off
    xlabel('Frequencies (Hz)');
    ylabel(sprintf('Coherence with %s', FeatureName));
    
    
    if ~isempty(FirstNonSigCoherenceFreq)
        if ~isempty(CoherencePeaksF)
            maxFreqInd = max(find(Freqs == max(CoherencePeaksF)), find(FirstNonSigCoherenceFreq==Freqs)) +2;
        else
            maxFreqInd = find(FirstNonSigCoherenceFreq==Freqs) +2;
        end
        subplot(1,3,3)
        shadedErrorBar(Freqs(1:maxFreqInd), Coherence_value(1:maxFreqInd),[(Coherence_up(1:maxFreqInd)-Coherence_value(1:maxFreqInd))'; (-Coherence_low(1:maxFreqInd)+Coherence_value(1:maxFreqInd))'], {'LineWidth',2,'Color','k'});
        hold on
        plot([Freqs(1) Freqs(maxFreqInd)], [0 0], 'r--', 'LineWidth',2);
        hold on
        plot(CoherencePeaksF, CoherencePeaks, 'go', 'MarkerSize',6,'MarkerFaceColor','g');
        hold on
        plot(MaxCoherence2, MaxCoherence1, 'ro', 'MarkerSize',6,'MarkerFaceColor','r');
        hold off
        xlabel('Frequencies (Hz)');
        ylabel(sprintf('Coherence with %s', FeatureName));
    end
end

Coherence.CoherencyT = CoherencyT;
Coherence.Coherence_value = Coherence_value;
Coherence.Coherence_up = Coherence_up;
Coherence.Coherence_low = Coherence_low;
Coherence.FirstNonSigCoherenceFreq = FirstNonSigCoherenceFreq;
Coherence.Info = Info;
Coherence.Info_up = Info_up;
Coherence.Info_low = Info_low;
Coherence.CoherencyT_filt =CoherencyT_filt;
Coherence.CoherencyT_xTimeDelay = CoherencyT_xTimeDelay;
Coherence.CoherencyT_DelayAtzero = CoherencyT_DelayAtzero;
Coherence.CoherencyT_WidthAtMaxPeak = CoherencyT_WidthAtMaxPeak;
Coherence.MaxCoherence1 = MaxCoherence1;
Coherence.MaxCoherence2 = MaxCoherence2;
Coherence.CoherencePeaks = CoherencePeaks;
Coherence.CoherencePeaksF = CoherencePeaksF;
Coherence.SecondCoherenceFreqCutOff = SecondCoherenceFreqCutOff;
Coherence.FeatureName = FeatureName;

end