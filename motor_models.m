addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'));
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'));
addpath(genpath('/Users/elie/Documents/CODE/LMC'));
DatFig=0; %Set to 1 to see input data figures for each cell
OutFig = 1;%Set to 1 to see output data figures for each cell

%% Listing datacells
%Filename = '59834_20190611_SSS_1-97.mat';
% Filename = '59834_20190610_SSS_1-130.mat';
Path = '/Users/elie/Documents/LMCResults/';
% Path = '/Users/elie/Documents/ManipBats/LMC/ResultsFiles/';
% Path = '/Users/elie/Google Drive/BatmanData/';

AllFiles = dir(fullfile(Path,'59834*.mat'));
Files2run = zeros(length(AllFiles),1);
for ff=1:length(AllFiles)
    if length(strfind(AllFiles(ff).name, '_'))==3
        Files2run(ff) = 1;
    end
end
CellsPath = AllFiles(logical(Files2run));

%% Running through cells to find the optimal time resolution of the neural response for acoustic feature predicion from the neural response
% here we calculate the coherency between the neural response and the
% acoustic features
TR=2; % 2ms is chosen as the Time resolution for the neural data
Fs = 1/(TR*10^-3); % the data is then sampled at the optimal frequency given the neural time resolution choosen
% find the closest power of 2 for the number of FFT window points that
% correspond to the Nyquist limit
Nyquist = Fs * 0.5;
nFFT = 2^ceil(log2(Nyquist));

NCells = length(CellsPath);
%Delay = nFFT/(2*Fs)*10^3;
 Delay=200;
% Lags = -Delay:Delay;
% Freqs = (0:ceil(length(Lags)/2)).* (2*Nyquist/length(Lags)); % Lags is a uneven number so F(i) = i*2*Nyquist/length(Lags)
MinCoherence4peaks = 0.01;
CoherencyT = cell(NCells,1);
CoherencyT_filt = cell(NCells,1);
CoherencyT_xTimeDelay = cell(NCells,1);
CoherencyT_DelayAtzero = nan(NCells,1);
CoherencyT_WidthAtMaxPeak = nan(NCells,1);
Freqs = cell(NCells,1);
% CoherencyF = cell(NCells,1);
Coherence = cell(NCells,1);
Coherence_low = cell(NCells,1);
Coherence_up = cell(NCells,1);
MaxCoherence = nan(NCells,2);
CoherencePeaks = cell(NCells,1);
CoherencePeaksF = cell(NCells,1);
FirstNonSigCoherenceFreq = nan(NCells,1);
SecondCoherenceFreqCutOff = nan(NCells,1);
Info = nan(NCells,1);
Info_low = nan(NCells,1);
Info_up = nan(NCells,1);
CellWithDurationIssue = [];

for cc=1:NCells % parfor
    fprintf(1, 'Cell %d/%d\n', cc, NCells)
    % load data
    Cell = load(fullfile(CellsPath(cc).folder,CellsPath(cc).name));
    
    % Number of vocalizations in the dataset
    if ~isfield(Cell, 'What')
        fprintf(1,'*** . Problem with Cell %d, no what field!! ****\n', cc)
        continue
    else
        IndVoc = find(contains(Cell.What, 'Voc'));
        NStims = length(IndVoc);
        StimDura = nan(NStims,1);
        for ss=1:NStims
            StimDura(ss) = round(length(Cell.BioSound{IndVoc(ss),2}.sound) ./(Cell.BioSound{IndVoc(ss),2}.samprate)*10^3);
        end
        if any(StimDura ~= Cell.Duration(IndVoc))
            CellWithDurationIssue = [CellWithDurationIssue cc];
            fprintf(1,'*** . Problem with Cell %d, duration inconcistency!! ****\n', cc)
            continue
        end
    end
        
    
    
    % Compute neural vectors
    % Neural Data loop
    % neural response is a vector that compile all spike counts starting
    % at -200ms (Delay) before stim onset and stop at 200ms (Delay) after
    % stim offset
    YPerStim = get_y_4Coherence(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),Delay,TR);
    Y = [YPerStim{:}]';    
    if sum(Y)<=2
        fprintf('Non Spiking cell, No calculation!!\n')
        continue
    end
    % Calculate acoustic features input to the models
    % acoustic data is a vector of the value of the acoustic feature sampled
    % at 1000Hz starting -200ms (Delay) before stim onset and stop at 200ms (Delay) after
    % stim offset
    DefaultVal = 0;%zero should be the default value for the amplitude, we know here that there is no sound
    XAmpPerStim = get_x_4coherence(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Delay,TR,DefaultVal,'amp');
    XAmp = [XAmpPerStim{:}]';
    
    % Calculate coherence and coherency between the signals
    [CoherencyT_unshifted, Freqs{cc}, CSR, CSR_up, CSR_low, stP] = multitapercoherence_JN([Y XAmp],nFFT,Fs);
    CoherencyT{cc} = fftshift(CoherencyT_unshifted);
    % Calculate the information on coherence
    %% normalize coherencies
    cStruct = struct;
    Freqs4Info = Freqs{cc};
    Coherence{cc} = CSR.^2;
    Coherence4Info = Coherence{cc};
    Coherence_up{cc} = CSR_up.^2;
    CoherenceUp4Info = Coherence_up{cc};
    
    closgn = sign(real(CSR_low));
    Coherence_low{cc} = (CSR_low.^2) .* closgn; %Coherence_low{cc} can be negative, multiply by sign after squaring
    CoherenceLow4Info = Coherence_low{cc};
    
    %% restrict frequencies analyzed to the requested cutoff and minimum frequency given the window size
    freqCutoff = Nyquist; % keep all frequencies up to Nyquist as of now
    if freqCutoff ~= -1
        findx = find(Freqs4Info < freqCutoff);
        eindx = max(findx);
        indx = 1:eindx;
        
        Freqs4Info = Freqs4Info(indx);
        Coherence4Info = Coherence4Info(indx);
        CoherenceUp4Info = CoherenceUp4Info(indx);
        CoherenceLow4Info = CoherenceLow4Info(indx);
    end 
    
    minFreq = Fs/nFFT;
    if minFreq > 0        
        findx = find(Freqs4Info >= minFreq);
        sindx = min(findx);
        Freqs4Info = Freqs4Info(sindx:end);
        Coherence4Info = Coherence4Info(sindx:end);
        CoherenceUp4Info = CoherenceUp4Info(sindx:end);
        CoherenceLow4Info = CoherenceLow4Info(sindx:end);        
    end
    
    %% if the clower goes below zero set all values to zero and keep track fo that first value of frequency for whoch the coherence is non significant
    cutoffIndex = find(CoherenceLow4Info < 0, 1, 'first');
    if (isempty(cutoffIndex))
        cutoffIndex = length(CoherenceLow4Info);        
    end
    freqCutoff = Freqs4Info(cutoffIndex);
    FirstNonSigCoherenceFreq(cc) = freqCutoff;
    
    %% compute information by integrating log of 1 - coherence        
    df = Freqs4Info(2) - Freqs4Info(1);
    cStruct.minFreq = minFreq;
    cStruct.freqCutoff = freqCutoff;
    Info(cc) = -df*sum(log2(1 - Coherence4Info(1:cutoffIndex)));
    Info_up(cc) = -df*sum(log2(1 - CoherenceUp4Info(1:cutoffIndex)));
    Info_low(cc) = -df*sum(log2(1 - CoherenceLow4Info(1:cutoffIndex)));
    
    %% Lowpass the coherency according to threshold FirstNonSigCoherenceFreq and find the width of max peak
    [z,p,k] = butter(6,FirstNonSigCoherenceFreq(cc)/(Fs/2),'low');
    sos_low = zp2sos(z,p,k);
    CoherencyT_filt{cc}=filtfilt(sos_low,1,CoherencyT{cc});
    CoherencyT_xTimeDelay{cc} = -(((nFFT/2)/Fs)*10^3):TR:(((nFFT/2-1))/Fs)*10^3; % Corresponding values in ms of the Delay for each value of CoherencyT
    [P,Locs] = findpeaks(CoherencyT_filt{cc});
    if ~isempty(Locs)
        [P,IndM] = max(P);
        Locs = Locs(IndM);
        CoherencyT_DelayAtzero(cc) = CoherencyT_xTimeDelay{cc}(Locs);
        CrossZero1 = CoherencyT_xTimeDelay{cc}(find(CoherencyT_filt{cc}(1:Locs)<=mean(CoherencyT_filt{cc}), 1, 'last')+1);
        CrossZero2 = CoherencyT_xTimeDelay{cc}(Locs + find(CoherencyT_filt{cc}(Locs+1:end)<=mean(CoherencyT_filt{cc}), 1, 'first')-1);
        CoherencyT_WidthAtMaxPeak(cc) = CrossZero2 - CrossZero1;
    end
    
    %% Find if there are other peaks in coherence than the initial one
    MaxCoherence(cc,1) = max(Coherence{cc});
    MaxCoherence(cc,2) = Freqs{cc}(Coherence{cc}==MaxCoherence(cc,1));
    [Coherence2ndPeaks, LocsC] = findpeaks(Coherence{cc}, 'MinPeakHeight', MinCoherence4peaks, 'MinPeakProminence',MinCoherence4peaks);
    LocsC(Coherence2ndPeaks==MaxCoherence(cc,1))=[];% eliminate the peak that corresponds to max value
    Coherence2ndPeaks(Coherence2ndPeaks==MaxCoherence(cc,1))=[]; % eliminate the peak that corresponds to max value
    Coherence2ndPeaks(Coherence_low{cc}(LocsC)<0)=[]; % eliminate peaks that are non-significant
    LocsC(Coherence_low{cc}(LocsC)<0)=[];% eliminate peaks that are non-significant
    CoherencePeaks{cc} = Coherence2ndPeaks';
    CoherencePeaksF{cc} = Freqs{cc}(LocsC)';
    SecondCoherenceFreqCutOff(cc) = max(CoherencePeaksF{cc});
    
    
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
%     CoherencyF{cc} = (FFT_XcorrAmpY)./(abs(FFT_AcorrAmp) .* abs(FFT_AcorrY)).^.5; % here we take abs(FFT) for autocorrelation to alleviate the effect of the phase (pick of the autocorrelation should be at zero phase and it's not when calculating the fft)
%     Coherence{cc} = (abs(CoherencyF{cc})).^2;
%     
%     % revert to time domain to find the coherency
%     CoherencyT{cc} = ifft(CoherencyF{cc});
    
    % Plot the value of coherency as a function of delay
    figure(3)
    clf
    subplot(1,3,1)
    plot(CoherencyT_xTimeDelay{cc}, CoherencyT_filt{cc}, 'LineWidth',2)
    if ~isempty(P)
        hold on
        plot(CoherencyT_DelayAtzero(cc),P,'ro', 'MarkerSize',16, 'MarkerFaceColor','r')
        text(max(CoherencyT_xTimeDelay{cc})/5, 4.5/5*P, sprintf('Delay at zero = %.2f ms', CoherencyT_DelayAtzero(cc)))
        text(max(CoherencyT_xTimeDelay{cc})/5, 4/5*P, sprintf('Width at y=mean(coherencyT) = %.2f ms', CoherencyT_WidthAtMaxPeak(cc)))
        hold off
    end
    xlabel('Time Delay')
    ylabel('Coherency')
%     subplot(1,2,2)
%     plot(Freqs, Coherence{cc}, '-k','LineWidth',2)
%     hold on
%     plot(Freqs, Coherence_up{cc}, '--r','LineWidth',2)
%     hold on
%     plot(Freqs, Coherence_low{cc}, '--b','LineWidth',2)
%     hold off
    subplot(1,3,2)
    shadedErrorBar(Freqs{cc}, Coherence{cc},[(Coherence_up{cc}-Coherence{cc})'; (-Coherence_low{cc}+Coherence{cc})'], {'LineWidth',2,'Color','k'})
    hold on
    plot([Freqs{cc}(1) Freqs{cc}(end)], [0 0], 'r--', 'LineWidth',2)
    hold on
    plot(CoherencePeaksF{cc}, CoherencePeaks{cc}, 'go', 'MarkerSize',6,'MarkerFaceColor','g')
    hold on
    plot(MaxCoherence(cc,2), MaxCoherence(cc,1), 'ro', 'MarkerSize',6,'MarkerFaceColor','r')
    MaxY = 0.4;
    ylim([-0.1 MaxY])
    text(Freqs{cc}(end)/3,0.8*MaxY,sprintf('Info = %.2f   InfoUp = %.2f    InfoLow = %.2f', Info(cc), Info_up(cc), Info_low(cc)))
    hold on
    text(Freqs{cc}(end)/3,0.9*MaxY,sprintf('MaxCoherence = %.2f', MaxCoherence(cc)))
    hold on
    text(Freqs{cc}(end)/3,0.85*MaxY,sprintf('FreqCutOff = %.2f Hz -> TimeResolution = %.2f ms', FirstNonSigCoherenceFreq(cc),1/FirstNonSigCoherenceFreq(cc)*10^3))
    hold off
    xlabel('Frequencies (Hz)')
    ylabel('Coherence')
    
    
    if ~isempty(FirstNonSigCoherenceFreq(cc))
        if ~isempty(CoherencePeaksF{cc})
            maxFreqInd = max(find(Freqs{cc} == max(CoherencePeaksF{cc})), find(FirstNonSigCoherenceFreq(cc)==Freqs{cc})) +2;
        else
            maxFreqInd = find(FirstNonSigCoherenceFreq(cc)==Freqs{cc}) +2;
        end
        subplot(1,3,3)
        shadedErrorBar(Freqs{cc}(1:maxFreqInd), Coherence{cc}(1:maxFreqInd),[(Coherence_up{cc}(1:maxFreqInd)-Coherence{cc}(1:maxFreqInd))'; (-Coherence_low{cc}(1:maxFreqInd)+Coherence{cc}(1:maxFreqInd))'], {'LineWidth',2,'Color','k'})
        hold on
        plot([Freqs{cc}(1) Freqs{cc}(maxFreqInd)], [0 0], 'r--', 'LineWidth',2)
        hold on
        plot(CoherencePeaksF{cc}, CoherencePeaks{cc}, 'go', 'MarkerSize',6,'MarkerFaceColor','g')
        hold on
        plot(MaxCoherence(cc,2), MaxCoherence(cc,1), 'ro', 'MarkerSize',6,'MarkerFaceColor','r')
        hold off
        xlabel('Frequencies (Hz)')
        ylabel('Coherence')
    end
%     keyboard
%     if Info(cc)>2
%         keyboard
%     end
end
    
save(fullfile(Path,'MotorModelsCoherency.mat'),'Delay','CoherencyT','CoherencyT_filt','CoherencyT_xTimeDelay','CoherencyT_DelayAtzero','CoherencyT_WidthAtMaxPeak','Freqs','Coherence','Coherence_low','Coherence_up','MaxCoherence','CoherencePeaks','CoherencePeaksF','FirstNonSigCoherenceFreq','SecondCoherenceFreqCutOff','Info','Info_low','Info_up','CellWithDurationIssue','CellsPath','TR','MinCoherence4peaks');


%% Plots results of coherence calculations for the population
% load(fullfile(Path,'MotorModelsCoherency.mat'))
figure(1)
clf
subplot(1,3,1)
histogram(MaxCoherence(:,1),'BinWidth',0.005,'FaceColor','k')
xlabel('Max coherence with sound amplitude')
ylabel('Number of cells')
hold on
v=vline(quantile(MaxCoherence(:,1), 0.25),'b--');
v.LineWidth = 2;
hold on
v=vline(quantile(MaxCoherence(:,1), 0.5),'g--');
v.LineWidth = 2;
hold on
v=vline(quantile(MaxCoherence(:,1), 0.75),'b--');
v.LineWidth = 2;

subplot(1,3,2)
plot(Info,MaxCoherence(:,1), 'o','Color','k','MarkerSize',6,'MarkerFaceColor','k')
xlabel('Information on Coherence with sound Amplitude(bits)')
ylabel('Max coherence with sound amplitude')

subplot(1,3,3)
histogram(Info,'BinWidth',0.05,'FaceColor','k')
xlabel('Information on Coherence with sound Amplitude(bits)')
ylabel('Number of Cells')
hold on
v=vline(quantile(Info, 0.25),'b--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.5),'g--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.75),'b--');
v.LineWidth = 2;

figure(2)
clf
subplot(2,5,1)
plot(Info, CoherencyT_DelayAtzero, 'o','Color','k','MarkerSize',6,'MarkerFaceColor','k')
xlabel('Information on Coherence with sound Amplitude(bits)')
ylabel('Phase of coherency (ms)')
hold on
v=vline(quantile(Info, 0.25),'b--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.5),'g--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.75),'b--');
v.LineWidth = 2;
hold on
h=hline(0,'r--');
h.LineWidth = 2;

subplot(2,5,2)
plot(Info, CoherencyT_WidthAtMaxPeak, 'o','Color','k','MarkerSize',6,'MarkerFaceColor','k')
xlabel('Information on Coherence with sound Amplitude(bits)')
ylabel('Time resolution of Coherency (ms)')
hold on
v=vline(quantile(Info, 0.25),'b--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.5),'g--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.75),'b--');
v.LineWidth = 2;

subplot(2,5,3)
plot(Info, FirstNonSigCoherenceFreq, 'o','Color','k','MarkerSize',6,'MarkerFaceColor','k')
xlabel('Information on Coherence with sound Amplitude(bits)')
ylabel('Max significant Frequency (Hz)')
hold on
v=vline(quantile(Info, 0.25),'b--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.5),'g--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.75),'b--');
v.LineWidth = 2;

subplot(2,5,4)
plot(Info, SecondCoherenceFreqCutOff, 'o','Color','k','MarkerSize',6,'MarkerFaceColor','k')
xlabel('Information on Coherence with sound Amplitude(bits)')
ylabel('Max 2nd peaks significant Frequency (Hz)')
hold on
v=vline(quantile(Info, 0.25),'b--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.5),'g--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.75),'b--');
v.LineWidth = 2;


subplot(2,5,5)
plot(Info, MaxCoherence(:,2), 'o','Color','k','MarkerSize',6,'MarkerFaceColor','k')
xlabel('Information on Coherence with sound Amplitude(bits)')
ylabel('Frequency of Max Coherence (Hz)')
hold on
v=vline(quantile(Info, 0.25),'b--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.5),'g--');
v.LineWidth = 2;
hold on
v=vline(quantile(Info, 0.75),'b--');
v.LineWidth = 2;
% Some points with very low values of info have high values of frequency of
% Max coherence, keep the plot focused on the majority of points
ylim([0 50])

subplot(2,5,6)
histogram(CoherencyT_DelayAtzero,'BinWidth',TR/2,'FaceColor','k')
xlabel('Phase of coherency (ms)')
ylabel('Number of cells')
hold on
v=vline(0, 'r:');
v.LineWidth = 2;

subplot(2,5,7)
histogram(CoherencyT_WidthAtMaxPeak,'BinWidth',TR/2,'FaceColor','k')
xlabel('Time resolution of Coherency (ms)')
ylabel('Number of cells')


subplot(2,5,8)
histogram(FirstNonSigCoherenceFreq, 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
xlabel('Max significant Frequency (Hz)')
ylabel('Number of Cells')

subplot(2,5,9)
histogram(SecondCoherenceFreqCutOff, 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
xlabel('Max 2nd peaks significant Frequency (Hz)')
ylabel('Number of Cells')

subplot(2,5,10)
histogram(MaxCoherence(:,2), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
% Some points with very low values of info have high values of frequency of
% Max coherence, keep the plot focused on the majority of points
xlim([0 50])
xlabel('Frequency of Max Coherence (Hz)')
ylabel('Number of Cells')

fprintf(1,'Cell with highest Info Value: %s', CellsPath(Info == max(Info)).name)

% Plot the values of the secondary peaks found for some cells
figure(4)
scatter([CoherencePeaksF{:}], [CoherencePeaks{:}],10, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor','k')
xlabel('Frequency of secondary peaks in Coherence')
ylabel('Values of Coherence')



% Restrict the dataset to Cells with Info >1 bit
GoodInfo = Info>=1;
% plot again the histograms for this subset of cells
figure(5)
clf
subplot(1,5,1)
histogram(CoherencyT_DelayAtzero(GoodInfo),'BinWidth',TR/2,'FaceColor','k')
xlabel('Phase of coherency (ms)')
ylabel('Number of cells with Info > 1bit')
hold on
v=vline(0, 'r:');
v.LineWidth = 2;

subplot(1,5,2)
histogram(CoherencyT_WidthAtMaxPeak(GoodInfo),'BinWidth',TR/2,'FaceColor','k')
xlabel('Time resolution of Coherency (ms)')
ylabel('Number of cells with Info > 1bit')


subplot(1,5,3)
histogram(FirstNonSigCoherenceFreq(GoodInfo), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
xlabel('Max significant Frequency (Hz)')
ylabel('Number of cells with Info > 1bit')

subplot(1,5,4)
histogram(SecondCoherenceFreqCutOff(GoodInfo), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
xlabel('Max 2nd peaks significant Frequency (Hz)')
ylabel('Number of Cells')

subplot(1,5,5)
histogram(MaxCoherence(find(GoodInfo),2), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
% Some points with very low values of info have high values of frequency of
% Max coherence, keep the plot focused on the majority of points
xlim([0 50])
xlabel('Frequency of Max Coherence (Hz)')
ylabel('Number of cells with Info > 1bit')


%% Motor Models parameters
% Assumption of stationarity over time

% Define the time resolution at which neural density estimates should be calculated
TRs = [75 100 125 150 175 200]; % time in ms


Win =300; %value in ms for the duration of the snippet of sound which
% acoustic features are used to predict the neural response at a given time
% t. This will be the size of the xaxis of the MRF.

Delay = 100;% ms
%The neural activity at time t is predicted by a time varying
% acoustic feature that starts at t-Delay ms

% parameters of the Poisson GLM Models
ParamModel.LAMBDARATIO=1e-4;
ParamModel.NUMLAMBDA=10;%25?
ParamModel.LINK='log';
ParamModel.DISTR='poisson';
% Determine a list of alpha (parameter that range the regularization
% betweeen ridge (L2, alpha=0) and lasso (L1, alpha =1))
ParamModel.Alpha=0.001; % STRFs are easier to interpret using ridge than using Lasso and deviances are similar.

%% Running through cells to find the optimal time resolution of the neural response for acoustic feature predicion from the neural response
% here we use a ridge regularization with an MSE optimization on the
% prediction of spike rate
NCells = sum(GoodInfo);
MSE_TR_Amp = cell(NCells,1);
MSE_TR_SpecMean = cell(NCells,1);
MSE_TR_Sal = cell(NCells,1);
Ypredict_Amp = cell(NCells,1);
Ypredict_SpecMean = cell(NCells,1);
Ypredict_Sal = cell(NCells,1);
Yval = cell(NCells,1);
MeanYTrain = cell(NCells,1);
TicToc = cell(NCells,1);
% load(fullfile(Path,'MotorModelsRidge.mat'));

for cc=1:NCells % parfor
    if ~isempty(MSE_TR_Amp{cc}) && ~isnan(MSE_TR_Amp{cc}(end)) && ~isempty(Yval{cc}) % This cell was already calculated
        fprintf(1, 'Cell %d/%d Already calculated\n',cc,NCells)
        continue
    end 
    Cell = load(fullfile(CellsPath(GoodInfo(cc)).folder,CellsPath(GoodInfo(cc)).name));
    
    
    % Number of vocalizations in the dataset
    if ~isfield(Cell, 'What')
        fprintf(1,'*** . Problem with Cell %s, no what field!! ****\n', CellsPath(GoodInfo(cc)).name)
        continue
    end
    IndVoc = find(contains(Cell.What, 'Voc') .* (Cell.Duration>20)); % No saliency calculated when duration <20ms
    NStims = length(IndVoc);
    
    %%  Split the set in a training 75% and testing set 25%
    AllRand = randperm(NStims);
    TrainSet = AllRand(1:round(NStims*0.75));
    ValSet = setdiff(AllRand, TrainSet);
    
    %% Calculate acoustic features input to the models
    % organize acoustic data as a matrix where each
    % column corresponds to amp env at t-100:Win:t+100, t being
    % time of neural window;
    % Acoustic feature loop, first column of biosound is microphone second is
    % piezo
    XSpecMeanPerStim = get_x(Cell.BioSound(IndVoc,1), Cell.Duration(IndVoc), Win, Delay,'SpectralMean');
    XSpecMeanTrain = [XSpecMeanPerStim{TrainSet}]';
    XSpecMeanTrain(isnan(XSpecMeanTrain))=nanmean(reshape(XSpecMeanTrain,numel(XSpecMeanTrain),1)); % Change the padding with NaN to the mean for the SpecMean for silence.
    XSpecMeanVal = [XSpecMeanPerStim{ValSet}]';
    XSpecMeanVal(isnan(XSpecMeanVal))=nanmean(reshape(XSpecMeanVal,numel(XSpecMeanVal),1)); % Change the padding with NaN to the mean for the SpecMean for silence.
    XSaliencyPerStim = get_x(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Win, Delay,'sal');
    XSaliencyTrain = [XSaliencyPerStim{TrainSet}]';
    XSaliencyTrain(isnan(XSaliencyTrain))=0; % Change the padding with NaN to zeros for the Saliency, we know Silence is not harmonic
    XSaliencyVal = [XSaliencyPerStim{ValSet}]';
    XSaliencyVal(isnan(XSaliencyVal))=0; % Change the padding with NaN to zeros for the Saliency, we know Silence is not harmonic
    XAmpPerStim = get_x(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Win, Delay,'amp');
    XAmpTrain = [XAmpPerStim{TrainSet}]';
    XAmpTrain(isnan(XAmpTrain))=0; % Change the padding with NaN to zeros for the amplitude, we know here that there is no sound
    XAmpVal = [XAmpPerStim{ValSet}]';
    XAmpVal(isnan(XAmpVal))=0; % Change the padding with NaN to zeros for the amplitude, we know here that there is no sound
    %         XSpecMedPerStim = get_x(BioSound(IndVoc,1), Duration(IndVoc), Win, TR, Delay,'Q2t');
    %         XSpecMed = [XSpecMedPerStim{:}]';
    
    
    
    %% Variable organization and running models
    MSE_TR_Amp{cc} = nan(1,length(TRs));
    MSE_TR_SpecMean{cc} = nan(1,length(TRs));
    MSE_TR_Sal{cc} =  nan(1,length(TRs));
    MeanYTrain{cc} = nan(1,length(TRs));
    Ypredict_Amp{cc} = cell(1,length(TRs));
    Ypredict_SpecMean{cc} = cell(1,length(TRs));
    Ypredict_Sal{cc} =  cell(1,length(TRs));
    Yval{cc} =  cell(1,length(TRs));
    TicToc{cc} = nan(1,length(TRs));
    Tr1=1;
    
    %%
    for tr = Tr1:length(TRs)
        TR = TRs(tr);
        fprintf(1,'Cell %d/%d Ridge models with Time resolution %d ms (%d/%d)\n', cc,NCells,TR, tr, length(TRs));
        %% Compute neural data vector
        TimerLoop =tic;
        % Neural Data loop
        % neural response is a vector that compile all spike counts starting
        % at -100ms (Delay) before stim onset and stop at 100ms (Delay) after
        % stim offset
        YPerStim = get_y(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),Win,Delay,TR);
        Y = [YPerStim{TrainSet}]';
        
        % save the values of Y with the smallest window for later model
        % evaluation at that smallest time resolution
        if tr==1
            Yval{cc}=[YPerStim{ValSet}]';
        end
        
        %% Plot of features and neuronal response if requested (BioSound,YPerStim,XPerStim, TR,Delay,F_high,FeatureName)
        if DatFig
            % Check the values of Y, to do calculations faster, we're going to
            % use a ridge regression, transforming the data with log to get
            % somewhat normal distributions, we don't want 0 values
            figure(2)
            clf
            subplot(1,2,1)
            histogram(Y)
            xlabel('Spike rate histogram Y')
            ylabel('# bins')
            subplot(1,2,2)
            histogram(log(Y))
            xlabel('Log Spike rate histogram log(Y)')
            ylabel('# bins')
            
            plotxyfeatures(BioSound(IndVoc,1),YPerStim,XSpecMeanPerStim,TR,Win,Delay,Duration(IndVoc), 50000,'Spectral Mean')
            plotxyfeatures(BioSound(IndVoc,2),YPerStim,XSaliencyPerStim,TR,Win,Delay,Duration(IndVoc), 10000,'Saliency')
            plotxyfeatures(BioSound(IndVoc,2),YPerStim,XAmpPerStim,TR,Win,Delay,Duration(IndVoc), 10000,'Amp')
            %         plotxyfeatures(BioSound(IndVoc,1),YPerStim,XSpecMedPerStim,TR,Delay,Duration(IndVoc), 50000,'Spectral Median')
        end
        
        %% Calculate the model
        % get rid of Nans
        Nan_Ind=find(isnan(Y));
        if ~isempty(Nan_Ind)
            keyboard
            Y(Nan_Ind) = [];
            XSpecMeanTrain(Nan_Ind,:) =[];
            XSaliencyTrain(Nan_Ind,:) =[];
            XAmpTrain(Nan_Ind,:) =[];
        end
        
        MeanYTrain{cc}(tr) = mean(Y);
    
    
        %% Run ridge regression on log transform of the data
        % Amp predicting Y
        [MSE_TR_Amp{cc}(tr),Ypredict_Amp{cc}{tr} ] = find_optimalTR(XAmpTrain,Y,XAmpVal,Yval{cc});
        % spectral mean predicting Y
        [MSE_TR_SpecMean{cc}(tr), Ypredict_SpecMean{cc}{tr}] = find_optimalTR(XSpecMeanTrain,Y,XSpecMeanVal, Yval{cc});
        % Pitch saliency predicting Y
        [MSE_TR_Sal{cc}(tr), Ypredict_Sal{cc}{tr}] = find_optimalTR(XSaliencyTrain,Y,XSaliencyVal,Yval{cc});
        TicToc{cc}(tr) = toc(TimerLoop);
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d) => done in %.2f minutes \n', cc,NCells,TR, tr, length(TRs), TicToc{cc}(tr)/60);
    end
    if OutFig
        figure(3)
        clf
        plot(TRs, MSE_TR_Amp{cc}, 'Linewidth',2)
        hold on
        plot(TRs, MSE_TR_SpecMean{cc}, 'Linewidth',2)
        hold on
        plot(TRs, MSE_TR_Sal{cc}, 'Linewidth',2)
        xlabel('Time resolution in ms')
        ylabel('Mean Squared Error')
        hold off
        
        figure(4)
        clf
        subplot(2,3,1)
        MAP = colormap();
        c=linspace(1,256,length(TRs));
        for tr = 1:length(TRs)
            scatter(Yval{cc}, Ypredict_Amp{cc}{tr},30,MAP(round(c(tr)),:),'filled')
            hold on
        end
        xlabel('Observed rate')
        ylabel('Predicted rate')
        title('Amplitude Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs)
        hold off
        YL = ylim;
        XL = xlim;
        ylim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        xlim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        
        subplot(2,3,2)
        for tr = 1:length(TRs)
            scatter(Yval{cc}, Ypredict_SpecMean{cc}{tr},30,MAP(round(c(tr)),:),'filled')
            hold on
        end
        xlabel('Observed rate')
        ylabel('Predicted rate')
        title('SpecMean Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs)
        hold off
        YL = ylim;
        XL = xlim;
        ylim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        xlim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        
        subplot(2,3,3)
        for tr = 1:length(TRs)
            scatter(Yval{cc}, Ypredict_Sal{cc}{tr},30,MAP(round(c(tr)),:),'filled')
            hold on
        end
        xlabel('Observed rate')
        ylabel('Predicted rate')
        title('Saliency Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs)
        hold off
        YL = ylim;
        XL = xlim;
        ylim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        xlim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        
        
        MSE_TR_Amp_local = nan(1,length(TRs));
        MSE_TR_SpecMean_local = nan(1,length(TRs));
        MSE_TR_Sal_local = nan(1,length(TRs));
        subplot(2,3,4)
        MAP = colormap();
        c=linspace(1,256,length(TRs));
        for tr = 1:length(TRs)
            scatter(Yval{cc}, (Ypredict_Amp{cc}{tr}-Yval{cc}).^2,30,MAP(round(c(tr)),:),'filled')
            hold on
            MSE_TR_Amp_local(tr) = mean((Ypredict_Amp{cc}{tr}-Yval{cc}).^2);
            line(xlim,MSE_TR_Amp_local(tr)*ones(2,1),'Color',MAP(round(c(tr)),:), 'LineWidth',2)
            hold on
        end
        xlabel('Observed rate')
        ylabel('Error2 on rate')
        title('Amplitude Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs)
        hold off
        
        
        subplot(2,3,5)
        for tr = 1:length(TRs)
            scatter(Yval{cc}, (Ypredict_SpecMean{cc}{tr}-Yval{cc}).^2,30,MAP(round(c(tr)),:),'filled')
            hold on
            MSE_TR_SpecMean_local(tr) = mean((Ypredict_SpecMean{cc}{tr}-Yval{cc}).^2);
            line(xlim,MSE_TR_SpecMean_local(tr)*ones(2,1),'Color',MAP(round(c(tr)),:), 'LineWidth',2)
            hold on
        end
        xlabel('Observed rate')
        ylabel('Error2 on log rate')
        title('SpecMean Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs)
        hold off
        
        subplot(2,3,6)
        for tr = 1:length(TRs)
            scatter(Yval{cc}, (Ypredict_Sal{cc}{tr}-Yval{cc}).^2,30,MAP(round(c(tr)),:),'filled')
            hold on
            MSE_TR_Sal_local(tr) = mean((Ypredict_Sal{cc}{tr}-Yval{cc}).^2);
            line(xlim,mean((Ypredict_Sal{cc}{tr}-Yval{cc}).^2)*ones(2,1),'Color',MAP(round(c(tr)),:), 'LineWidth',2)
            hold on
        end
        xlabel('Observed rate')
        ylabel('Error2 on log rate')
        title('Saliency Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs)
        hold off
        
        figure(3)
        hold on
        plot(TRs, MSE_TR_Amp_local, 'Linewidth',2, 'LineStyle','--')
        hold on
        plot(TRs, MSE_TR_SpecMean_local, 'Linewidth',2,'LineStyle','--')
        hold on
        plot(TRs, MSE_TR_Sal_local, 'Linewidth',2,'LineStyle','--')
        legend({'Amp' 'SpectralMean' 'Saliency' 'AmpMe' 'SpectralMeanMe' 'SaliencyMe'})
        hold off
        
        figure(5)
        clf
        plot(TRs,MeanYTrain{cc}-mean(Yval{cc}))
        xlabel('Time resolution (ms)')
        ylabel('Mean rate difference Training-testing')
    end
%     keyboard
end
save(fullfile(Path,'MotorModelsRidge.mat'), 'NCells', 'MSE_TR_Amp', 'MSE_TR_SpecMean','MSE_TR_Sal','Ypredict_Amp','Ypredict_SpecMean','Ypredict_Sal','Yval','MeanYTrain','TicToc','CellsPath');


%% Plot the results of the time resolution optimization using ridge regression
load(fullfile(Path,'MotorModelsRidge.mat'));
% Plot the zscored average MSE accross cells
MSE_TR_Amp_zs = nan(NCells,length(TRs));
MSE_TR_SpecMean_zs = nan(NCells,length(TRs));
MSE_TR_Sal_zs = nan(NCells,length(TRs));
% also the MSE as a proportion of the average rate of the training correspoding portion
% of the validating set (same temporal resolution)
MSE0 = nan(NCells,length(TRs));
R2_TR_Amp = nan(NCells,length(TRs));
R2_TR_SpecMean = nan(NCells,length(TRs));
R2_TR_Sal = nan(NCells,length(TRs));
for cc=1:142
    if isempty(MSE_TR_Amp{cc})
        continue
    end
    MSE_TR_Amp_zs(cc,:) = zscore(MSE_TR_Amp{cc});
    MSE_TR_SpecMean_zs(cc,:) = zscore(MSE_TR_SpecMean{cc});
    MSE_TR_Sal_zs(cc,:) = zscore(MSE_TR_Sal{cc});
    for tr=1:length(TRs)
        MSE0(cc,tr) = mean((Yval{cc}-MeanYTrain{cc}(tr)).^2);
    end
    R2_TR_Amp(cc,:) = (MSE0(cc,:)-MSE_TR_Amp{cc})./MSE0(cc,:);
    R2_TR_SpecMean(cc,:) = (MSE0(cc,:)-MSE_TR_SpecMean{cc})./MSE0(cc,:);
    R2_TR_Sal(cc,:) = (MSE0(cc,:)-MSE_TR_Sal{cc})./MSE0(cc,:);
end


figure(7)
ColorCode = get(groot, 'DefaultAxesColorOrder');
subplot(1,2,1)
legend('AutoUpdate', 'on')
plot(TRs, nanmean(MSE_TR_Amp_zs), 'LineWidth',2, 'Color',ColorCode(1,:), 'DisplayName','Amplitude')
hold on
plot(TRs, nanmean(MSE_TR_SpecMean_zs), 'LineWidth',2, 'Color',ColorCode(2,:),'DisplayName','SpectralMean')
hold on
plot(TRs, nanmean(MSE_TR_Sal_zs), 'LineWidth',2, 'Color',ColorCode(3,:),'DisplayName','Saliency')
legend('AutoUpdate', 'off')
hold on
shadedErrorBar(TRs, nanmean(MSE_TR_Amp_zs),nanstd(MSE_TR_Amp_zs)./(sum(~isnan(MSE_TR_Amp_zs))).^0.5, {'Color',ColorCode(1,:)})
hold on
shadedErrorBar(TRs, nanmean(MSE_TR_SpecMean_zs),nanstd(MSE_TR_SpecMean_zs)./(sum(~isnan(MSE_TR_SpecMean_zs))).^0.5, {'Color', ColorCode(2,:)})
hold on
shadedErrorBar(TRs, nanmean(MSE_TR_Sal_zs),nanstd(MSE_TR_Sal_zs)./(sum(~isnan(MSE_TR_Sal_zs))).^0.5, {'Color', ColorCode(3,:)})
hold on
plot(TRs, nanmean(MSE_TR_Amp_zs), 'LineWidth',2, 'Color',ColorCode(1,:))
hold on
plot(TRs, nanmean(MSE_TR_SpecMean_zs), 'LineWidth',2, 'Color',ColorCode(2,:))
hold on
plot(TRs, nanmean(MSE_TR_Sal_zs), 'LineWidth',2, 'Color',ColorCode(3,:))
hold off
xlabel('Time resolution in ms')
ylabel('zscored Mean Squared Error')

subplot(1,2,2)
legend('AutoUpdate', 'on')
plot(TRs, nanmean(R2_TR_Amp), 'LineWidth',2, 'Color',ColorCode(1,:), 'DisplayName','Amplitude')
hold on
plot(TRs, nanmean(R2_TR_SpecMean), 'LineWidth',2, 'Color',ColorCode(2,:),'DisplayName','SpectralMean')
hold on
plot(TRs, nanmean(R2_TR_Sal), 'LineWidth',2, 'Color',ColorCode(3,:),'DisplayName','Saliency')
hold on
legend('AutoUpdate', 'off')
shadedErrorBar(TRs, nanmean(R2_TR_Amp),nanstd(R2_TR_Amp)./(sum(~isnan(R2_TR_Amp))).^0.5, {'Color',ColorCode(1,:)})
hold on
shadedErrorBar(TRs, nanmean(R2_TR_SpecMean),nanstd(R2_TR_SpecMean)./(sum(~isnan(R2_TR_SpecMean))).^0.5, {'Color', ColorCode(2,:)})
hold on
shadedErrorBar(TRs, nanmean(R2_TR_Sal),nanstd(R2_TR_Sal)./(sum(~isnan(R2_TR_Sal))).^0.5, {'Color', ColorCode(3,:)})
hold on

plot(TRs, nanmean(R2_TR_Amp), 'LineWidth',2, 'Color',ColorCode(1,:))
hold on
plot(TRs, nanmean(R2_TR_SpecMean), 'LineWidth',2, 'Color',ColorCode(2,:))
hold on
plot(TRs, nanmean(R2_TR_Sal), 'LineWidth',2, 'Color',ColorCode(3,:))
hold off
xlabel('Time resolution in ms')
ylabel('R2')
hold off

figure()
subplot(1,3,1)
imagesc(R2_TR_Amp)
set(gca,'XTickLabel', TRs)
xlabel('Time Resolution (ms)')
ylabel('R2 Amplitude')
colorbar()
caxis([0 0.05])

subplot(1,3,2)
imagesc(R2_TR_SpecMean)
set(gca,'XTickLabel', TRs)
xlabel('Time Resolution (ms)')
ylabel('R2 SpectralMean')
colorbar()
caxis([0 0.05])

subplot(1,3,3)
imagesc(R2_TR_Sal)
set(gca,'XTickLabel', TRs)
xlabel('Time Resolution (ms)')
ylabel('R2 Saliency')
colorbar()
caxis([0 0.05])



%%         %% Run ridge GLM Poisson on acoustic features

% spectral mean predicting Y
 [Mdl_SpecMean_local]=fitrlinear([XAmpTrain XSpecMeanTrain],log(Y),'Regularization', 'ridge', 'Learner','leastsquares', 'CrossVal','on');
        % find the model with the minimum of deviance (best lambda)
        [BestDevSpecMean(tr),BestModSpecMean] = min(FitInfo_SpecMean.Deviance);
        BSpecMean(tr,2:end) = BSpecMean_local(:,BestModSpecMean);
        BSpecMean(tr,1) = FitInfo_SpecMean.Intercept(BestModSpecMean);
        
        % pitch saliency predicting Y
        [BSal_local, FitInfo_Sal]=lassoglm([XAmpTrain XSaliencyTrain],Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',1,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevSal(tr),BestModSal] = min(FitInfo_Sal.Deviance);
        BSal(tr,2:end) = BSal_local(:,BestModSal);
        BSal(tr,1) = FitInfo_Sal.Intercept(BestModSal);
        
        
        %         % spectral median predicting Y
        %         [BSpecMed, FitInfo_SpecMed]=lassoglm(XSpecMed,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
        %         % find the model with the minimum of deviance (best lambda)
        %         [BestDevSpecMed(tr,dd),BestModSpecMed] = min(FitInfo_SpecMean.Deviance);
        
        % amplitude predicting Y, is a null model for the former 2 models
        [BAmp_local, FitInfo_Amp]=lassoglm(XAmpTrain,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',1,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevAmp(tr),BestModAmp] = min(FitInfo_Amp.Deviance);
        
        BAmp(tr,2:end) = BAmp_local(:,BestModAmp);
        BAmp(tr,1) = FitInfo_Amp.Intercept(BestModAmp);
        
        % null model
        MDL_null=fitglm(ones(size(XAmpTrain,1),1),Y,'Distribution',ParamModel.DISTR,'Link',ParamModel.LINK,'Intercept',false);
        % find the model with the minimum of deviance (best lambda)
        %     [BestDevNull(tr),BestModNull] = min(FitInfo_Amp.Deviance);
        %
        %     BNull(tr,2:end) = Bnull_local(:,BestModNull);
        %     BNull(tr,1) = FitInfo_null.Intercept(BestModNull);
        
        BestDevNull(tr) = MDL_null.Deviance;
        BNull(tr) = MDL_null.Coefficients.Estimate;

        if DatFig
            TimeBinsX = -Delay : (Win-Delay);
            figure()
            plot(TimeBinsX,BSpecMean(:,BestModSpecMean),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Spectral Mean Dev = %.1f',BestDevSpecMean(tr)))
            hold on
            plot(TimeBinsX,BSal(:,BestModSal),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Saliency Dev = %.1f',BestDevSal(tr)))
            hold on
            %         plot(TimeBinsX,BSpecMed(:,BestModSpecMed),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Spectral Median Dev = %.1f',BestDevSpecMed(tr,dd)))
            hold on
            plot(TimeBinsX,BAmp(:,BestModAmp),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Amp Dev = %.1f',BestDevAmp(tr)))
            
            
            legend('show')
            % XTick = get(gca, 'XTickLabel');
            % XTick = cellfun(@str2double, XTick) * Win;
            % set(gca,'XTickLabel',XTick)
            xlabel('Time (ms)')
            title(sprintf('Poisson Ridge regression on Acoustic Features Delay = %dms Time resolution=%dms', Delay, TR))
            hold off
            pause(1)
        end
        TicToc(tr) = toc(TimerLoop);
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d) => done in %.2f minutes \n', cc,NCells,TR, tr, length(TRs), TicToc(tr)/60);
   
    if DatFig
        
        figure()
        ColorCode = get(groot, 'DefaultAxesColorOrder');
        subplot(1,2,1)
        plot(TRs,BestDevSpecMean', 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        hold on
        plot(TRs,BestDevNull, 'Color',ColorCode(4,:),'LineWidth',2, 'LineStyle','--')
        legend({'Amp + SpecMean' 'Amp + Sal' 'Amp' 'Null'})
        legend('Autoupdate','off')
        ylabel('Deviance')
        xlabel('Time resolution in ms')
        YLim = get(gca,'YLim');
        
        subplot(1,2,2)
        plot(TRs,BestDevSpecMean', 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        hold on
        plot(TRs,BestDevNull, 'Color',ColorCode(4,:),'LineWidth',2, 'LineStyle','--')
        ylabel('Deviance')
        xlabel('Time resolution in ms')
        xlim([0 2])
        ylim([0.9*YLim(2) YLim(2)])
        hold off
        
        
        figure()
        subplot(1,2,1)
        plot(TRs, BestDevAmp-BestDevSpecMean , 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevAmp - BestDevSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevNull -BestDevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
        ylabel('Difference of Deviance')
        xlabel('Time resolution in ms')
        
        subplot(1,2,2)
        plot(TRs, (BestDevAmp-BestDevSpecMean)./BestDevAmp , 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs,(BestDevAmp - BestDevSal)./BestDevAmp,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs,(BestDevNull -BestDevAmp)./BestDevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
        ylabel('R2')
        xlabel('Time resolution in ms')
        
        % Plotting betas
        figure()
        subplot(2,3,1:3)
        bar(TRs,exp([BNull BAmp(:,1) BSpecMean(:,1) BSal(:,1)])*1000)
        xlabel('Time resolution in ms')
        ylabel('exp(intercept) Hz')
        legend({'Null' 'Amp' 'Amp + SpecMean' 'Amp + Sal'},'Location','southoutside','NumColumns',4)
        subplot(2,3,4)
        imagesc(BAmp(:,2:end))
        set(gca,'YTick',1:length(TRs),'YTickLabel',TRs)
        ylabel('Time resolution in ms')
        xlabel('Time alligned to predicted Y (ms)')
        set(gca,'XTick',1:50:Win,'XTickLabel',-Delay:50:(Win-Delay-1))
        title('Amplitude Model')
        colorbar()
        subplot(2,3,5)
        imagesc(BSpecMean(:,2:end))
        set(gca,'YTick',1:length(TRs),'YTickLabel',TRs)
        ylabel('Time resolution in ms')
        xlabel('Time alligned to predicted Y (ms)')
        set(gca,'XTick',1:100:2*Win,'XTickLabel',[-Delay:100:(Win-Delay-1) -Delay:100:(Win-Delay-1)])
        title('Amplitude + SpecMean')
        colorbar()
        subplot(2,3,6)
        imagesc(BSal(:,2:end))
        set(gca,'YTick',1:length(TRs),'YTickLabel',TRs)
        ylabel('Time resolution in ms')
        xlabel('Time alligned to predicted Y (ms)')
        set(gca,'XTick',1:100:2*Win,'XTickLabel',[-Delay:100:(Win-Delay-1) -Delay:100:(Win-Delay-1)])
        title('Amplitude + Saliency')
        colorbar()
    end
    
    
    MotorModels{cc}.TRs = TRs;
    MotorModels{cc}.DevAmp = BestDevAmp;
    MotorModels{cc}.DevAmpSpecMean = BestDevSpecMean;
    MotorModels{cc}.DevAmpSal = BestDevSal;
    MotorModels{cc}.DevNull = BestDevNull;
    MotorModels{cc}.BAmp = BAmp;
    MotorModels{cc}.BAmpSpecMean = BSpecMean;
    MotorModels{cc}.BAmpSal = BSal;
    MotorModels{cc}.BNull = BNull;
    MotorModels{cc}.TicToc = TicToc;
    MotorModels{cc}.MeanY = MeanY;
        

%%
save(fullfile(Path,'MotorModelsAllCells'),'MotorModels','CellsPath');
%% Plot issues with no convergence
load(fullfile(Path,'MotorModelsAllCells'),'MotorModels','CellsPath');
TicToc = nan(NCells*11,1);
MeanY = nan(NCells*11,1);
TRs = nan(NCells*11,1);
BestDevAmp = nan(NCells*11,1);
BestDevSpecMean = nan(NCells*11,1);
BestDevSal = nan(NCells*11,1);
BestDevNull = nan(NCells*11,1);

F10=figure(10);
F20 = figure(20);
ColorCode = get(groot, 'DefaultAxesColorOrder');
ii=0;
for cc=1:NCells
    if ~isempty(MotorModels{cc})
        for tt=1:length(MotorModels{cc}.TicToc)
            ii=ii+1;
            TicToc(ii) = MotorModels{cc}.TicToc(tt);
            MeanY(ii) = MotorModels{cc}.MeanY(tt);
            TRs(ii) = MotorModels{cc}.TRs(tt);
            BestDevAmp(ii) = MotorModels{cc}.DevAmp(tt);
            BestDevSpecMean(ii) = MotorModels{cc}.DevAmpSpecMean(tt);
            BestDevSal(ii) = MotorModels{cc}.DevAmpSal(tt);
            BestDevNull(ii) = MotorModels{cc}.DevNull(tt);
        end
        
        figure(10)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevAmpSpecMean', 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevAmpSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevNull, 'Color',ColorCode(4,:),'LineWidth',2, 'LineStyle','--')
        if cc==1
            legend({'Amp + SpecMean' 'Amp + Sal' 'Amp' 'Null'})
            legend('Autoupdate','off')
            ylabel('Deviance')
            xlabel('Time resolution in ms')
        end
        
        
        figure(20)
        subplot(1,2,1)
        hold on
        plot(MotorModels{cc}.TRs, MotorModels{cc}.DevAmp-MotorModels{cc}.DevAmpSpecMean , 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevAmp - MotorModels{cc}.DevAmpSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevNull -MotorModels{cc}.DevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        if cc==1
            legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
            legend('Autoupdate','off')
            ylabel('Difference of Deviance')
            xlabel('Time resolution in ms')
        end
        
        subplot(1,2,2)
        hold on
        plot(MotorModels{cc}.TRs, (MotorModels{cc}.DevAmp-MotorModels{cc}.DevAmpSpecMean)./MotorModels{cc}.DevAmp , 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,(MotorModels{cc}.DevAmp - MotorModels{cc}.DevAmpSal)./MotorModels{cc}.DevAmp,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,(MotorModels{cc}.DevNull -MotorModels{cc}.DevAmp)./MotorModels{cc}.DevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        if cc==1
            legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
            legend('Autoupdate','off')
            ylabel('R2')
            xlabel('Time resolution in ms')
        end
        
    end
end

figure()
subplot(1,2,1)
scatter(MeanY,TicToc/(60*60), 'MarkerFaceColor','k')
xlabel('MeanRate')
ylabel('Model Time consumption (h)')
subplot(1,2,2)
scatter(TRs,TicToc/(60*60), 'MarkerFaceColor','k')
xlabel('Models Time Resolution (ms)')
ylabel('Model Time consumption (h)')


figure()
ColorCode = get(groot, 'DefaultAxesColorOrder');
%subplot(1,2,1)
scatter(TRs,BestDevSpecMean, 'MarkerFaceColor',ColorCode(1,:), 'MarkerEdgeColor',ColorCode(1,:))
hold on
scatter(TRs,BestDevSal,'MarkerFaceColor',ColorCode(2,:), 'MarkerEdgeColor',ColorCode(2,:))
hold on
scatter(TRs,BestDevAmp, 'MarkerFaceColor',ColorCode(3,:), 'MarkerEdgeColor',ColorCode(3,:))
hold on
scatter(TRs,BestDevNull, 'MarkerFaceColor',ColorCode(4,:), 'MarkerEdgeColor',ColorCode(4,:))
legend({'Amp + SpecMean' 'Amp + Sal' 'Amp' 'Null'})
legend('Autoupdate','off')
ylabel('Deviance')
xlabel('Time resolution in ms')




figure()
subplot(1,2,1)
scatter(TRs, BestDevAmp-BestDevSpecMean , 'MarkerFaceColor',ColorCode(1,:), 'MarkerEdgeColor',ColorCode(1,:))
hold on
scatter(TRs,BestDevAmp - BestDevSal,'MarkerFaceColor',ColorCode(2,:), 'MarkerEdgeColor',ColorCode(2,:))
hold on
scatter(TRs,BestDevNull -BestDevAmp, 'MarkerFaceColor',ColorCode(3,:), 'MarkerEdgeColor',ColorCode(3,:))
legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
ylabel('Difference of Deviance')
xlabel('Time resolution in ms')

subplot(1,2,2)
scatter(TRs, (BestDevAmp-BestDevSpecMean)./BestDevAmp , 'MarkerFaceColor',ColorCode(1,:), 'MarkerEdgeColor',ColorCode(1,:))
hold on
scatter(TRs,(BestDevAmp - BestDevSal)./BestDevAmp,'MarkerFaceColor',ColorCode(2,:), 'MarkerEdgeColor',ColorCode(2,:))
hold on
scatter(TRs,(BestDevNull -BestDevAmp)./BestDevAmp,'MarkerFaceColor',ColorCode(3,:), 'MarkerEdgeColor',ColorCode(3,:))
legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
ylabel('R2')
xlabel('Time resolution in ms')
%%
        





% %% Old plots of deviances
% CLIM = nan(4,2);
% figure()
% imagesc(Delay,TRs,BestDevSal)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance pitch saliency model')
% colorbar()
% CLIM(2,:) = get(gca, 'clim');
% 
% % figure()
% % imagesc(Delays,TRs,BestDevSpecMed)
% % xlabel('Delays to voc onset in ms')
% % ylabel('Time resolution in ms')
% % title('Deviance spectral median model')
% 
% figure()
% imagesc(Delay,TRs,BestDevAmp)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance Amplitude model')
% colorbar()
% CLIM(3,:) = get(gca, 'clim');
% 
% figure()
% imagesc(Delay,TRs,BestDevNull)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance Null model')
% colorbar()
% CLIM(4,:) = get(gca, 'clim');
% 
% 
% figure()
% subplot(2,3,1)
% imagesc(BestDevSpecMean)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance spectral mean model')
% colorbar()
% %caxis([0 1000])
% caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% 
% subplot(2,3,2)
% imagesc(BestDevSal)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance pitch saliency model')
% colorbar()
% caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% %caxis([0 1000])
% 
% % figure()
% % imagesc(Delays,TRs,BestDevSpecMed)
% % xlabel('Delays to voc onset in ms')
% % ylabel('Time resolution in ms')
% % title('Deviance spectral median model')
% 
% subplot(2,3,3)
% imagesc(BestDevAmp)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance Amplitude model')
% colorbar()
% caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% %caxis([0 1000])
% 
% subplot(2,3,5)
% imagesc(BestDevNull)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance Null model')
% colorbar()
% caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% 
% figure()
% subplot(1,3,1)
% imagesc(BestDevAmp - BestDevSpecMean)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Goodness of fit: Deviance Amp - SpecMean')
% colorbar()
% % caxis([0 1000])
% % %caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% 
% subplot(1,3,2)
% imagesc(BestDevAmp - BestDevSal)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Goodness of fit: Deviance Amp - pitch saliency')
% colorbar()
% % %caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% % caxis([0 1000])
% 
% % figure()
% % imagesc(Delays,TRs,BestDevSpecMed)
% % xlabel('Delays to voc onset in ms')
% % ylabel('Time resolution in ms')
% % title('Deviance spectral median model')
% 
% subplot(1,3,3)
% imagesc(BestDevNull - BestDevAmp)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Goodness of fit: Deviance Null -  Amplitude')
% colorbar()
% %caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% caxis([0 1000])

% %% Run ridge GLM Poisson on change of power
% B = cell(length(Fhigh_powers),1);
% Dev_Env = nan(length(Fhigh_powers),1);
% FitInfo_Env = cell(length(Fhigh_powers),1);
% figure()
% for ff=1:length(Fhigh_powers)
%     X_local = abs(diff(X_Pow{ff}')');
% %     [B_Env{ff}, Dev_Env(ff)]=glmfit(X_Pow{ff},Y,ParamModel.DISTR,'link',ParamModel.LINK);
%     [B{ff}, FitInfo_Env{ff}]=lassoglm(X_local,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
%     % find the model with the minimum of deviance (best lambda)
%     [BestDevSal,BestModSal] = min(FitInfo_Env{ff}.Deviance);
%     plot(TimeBinsX(2:end),B{ff}(:,BestModSal),'Color',ColorCode(ff+1,:),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('FHigh %d Hz Dev = %.1f',Fhigh_powers(ff),BestDevSal))
%     hold on
% end
% legend('show')
% % XTick = get(gca, 'XTickLabel');
% % XTick = cellfun(@str2double, XTick) * Win;
% % set(gca,'XTickLabel',XTick)
% xlabel('Time (ms)')
% title('Poisson Ridge regression on change of Power')
%
%
% %% run logistic model
% Y01 = Y>0;
% B_Env01 = cell(length(Fhigh_powers),1);
% Dev_Env01 = nan(length(Fhigh_powers),1);
% FitInfo_Env01 = cell(length(Fhigh_powers),1);
% figure()
% for ff=1:length(Fhigh_powers)
% %     [B_Env01{ff}, Dev_Env01(ff)]=glmfit(X_Pow{ff},Y01,'binomial','link','logit');
%     [B_Env01{ff}, FitInfo_Env01{ff}]=lassoglm(X_Pow{ff},Y01,'binomial','Alpha', ParamModel.Alpha,'Link','logit','NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
%     % find the model with the minimum of deviance (best lambda)
%     [BestDevSal,BestModSal] = min(FitInfo_Env01{ff}.Deviance);
%     plot(TimeBinsX,B{ff}(:,BestModSal),'Color',ColorCode(ff+1,:),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('FHigh %d Hz Dev = %.1f',Fhigh_powers(ff),BestDevSal))
%     hold on
% end
% legend('show')
% XTick = get(gca, 'XTickLabel');
% XTick = cellfun(@str2double, XTick) * TR;
% set(gca,'XTickLabel',XTick)
% xlabel('Time (ms)')
% title('Logistic Ridge regression on Power')
% %     [B_Env{ff}, FitInfo_Env]=lassoglm(X,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
%
% % lassoPlot(B_Env,FitInfo_Env,'PlotType','Lambda','XScale','log')
%



%% INTERNAL FUNCTIONS

function [XPerStim] = get_x(BioSound, Duration, Win,Delay,Feature)
% calculate the cell array of acoustic data for each vocalization. For each
% voc, XPerStim is a matrix where each column correspond to the values of
% the acoustic feature in the window Win starting Delay ms before time t of
% the neural response Y. Same as the neural response, the window slides by
% a pace of 1ms from one column of the matrix to the next, starting
% Win-Delay before the vocalization onset.

% Matrix of the acoustic features
XPerStim = cell(1,length(Duration));
for stim = 1:length(Duration)
    % Time points of the corresponding neural response
    TimeBinsY = -(Win-Delay) : 1 : (Delay + Duration(stim));
    
    % Get ready the stim acoustic features that was sampled at 1000Hz
    FeatureVal = BioSound{stim}.(sprintf('%s',Feature));
    %         FeatureValResampled = resample(FeatureVal,1,TR); % Win in ms so resampling from 1 value per ms to 1/Win value per ms or 1 value per win ms
    %         PadFeature4Conv = nanmean(FeatureVal)*ones(1,(length(Expwav)-1)/2);
    %         FeatureValResampled = conv([PadFeature4Conv FeatureVal PadFeature4Conv], Expwav, 'valid');
    ZPaddFeature = nan(1,Win + Duration(stim) + Win);
    ZPaddFeature(Win+(1:length(FeatureVal)))=FeatureVal;
    
    XPerStim{stim} = zeros(Win,length(TimeBinsY)-1);
    for tt=1:(length(TimeBinsY)-1)
        XPerStim{stim}(:,tt) = ZPaddFeature((1:Win) + tt-1);
    end
end

end

function [XPerStim] = get_x_4coherence(BioSound, Duration, Delay,TR,DefaultVal,Feature)
% calculate the cell array of acoustic data for each cell. For each
% voc, XPerStim is a vector where each column correspond to the values of
% the acoustic feature in the window Win starting Delay ms before time onset
% of the vocalization.
Fs = 1/(TR.*10^-3);
% Vectors of the acoustic features
XPerStim = cell(1,length(Duration));
for stim = 1:length(Duration)
    if round(length(BioSound{stim}.sound)/BioSound{stim}.samprate*10^3) ~= Duration(stim)
        keyboard
    end
    % Get ready an output vector for the stim acoustic features that was sampled at 1000Hz
    XPerStim_temp = DefaultVal.*ones(1,round(1000 .* (Delay + Duration(stim) + Delay).*10^-3));
    FeatureVal = BioSound{stim}.(sprintf('%s',Feature));
    XPerStim_temp(Delay+(1:length(FeatureVal)))=FeatureVal;
    XPerStim{stim} = resample(XPerStim_temp, Fs, 1000); 
end

end

function [YPerStim] = get_y(SAT, Duration, Win,Delay,TR)
% Calculate the time varying rate applying a gaussian window TR on the
% spike pattern. The spike pattern considered starts Win-Delay ms
% before the onset of the vocalization and stops Delay ms after the
% offset of the vocalization
YPerStim = cell(1,length(Duration));
% Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
nStd =max(Duration) + Win-Delay + Delay; % before set as 4
Tau = (TR/2);
T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
Expwav = Expwav./sum(Expwav);

% Loop through the stimuli and fill in the matrix
for stim=1:length(Duration)
    % Time slots for the neural response
    TimeBinsY = -(Win-Delay) : (Delay + Duration(stim));
    SpikePattern = zeros(1,length(TimeBinsY)-1);
    for isp = 1:length(SAT{stim})
        SpikeInd = round(SAT{stim}(isp));
        if (SpikeInd>=-(Win-Delay)) && (SpikeInd<(Delay + Duration(stim)))
            SpikePattern(SpikeInd + (Win-Delay) +1) = SpikePattern(SpikeInd + (Win-Delay) +1) +1;
        end
    end
    YPerStim{stim} = conv(SpikePattern, Expwav,'same');
%     % change zero values for the smallest value under matlab.
%     if sum(YPerStim{stim}==0)
%         MinData = min(YPerStim{stim}(YPerStim{stim} ~=0));
%         if ~isempty(MinData)
%             YPerStim{stim}(YPerStim{stim}==0)=min(MinData,realmin('double'));
%         else
%             YPerStim{stim}(YPerStim{stim}==0)=realmin('double');
%         end
%     end
    
    % Make sure that the output sum(Y) = input sum(Y)
    if (~sum(YPerStim{stim})==0) && ~(sum(SpikePattern)==0)
        YPerStim{stim} = YPerStim{stim}./sum(YPerStim{stim}).*sum(SpikePattern);
    end
    
    
    %         % Time slots for the neural response
    %         TimeBinsY = -(Delay) : TR: (Delay + Duration(stim));
    %         YPerStim{stim} = nan(1,length(TimeBinsY)-1);
    %         for tt=1:(length(TimeBinsY)-1)
    %             % Find the number of spikes
    %             YPerStim{stim}(tt) = sum( (SAT{stim}>=TimeBinsY(tt)) .* (SAT{stim}<TimeBinsY(tt+1)));
    %         end
end

end

function [YPerStim] = get_y_4Coherence(SAT, Duration,Delay,TR)
% Calculate the time varying rate applying a gaussian window TR on the
% spike pattern. The spike pattern considered starts -Delay ms
% before the onset of the vocalization and stops Delay ms after the
% offset of the vocalization
YPerStim = cell(1,length(Duration));
% Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
nStd =max(Duration) + Delay + Delay; % before set as 4
Tau = (TR/2);
T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
Expwav = Expwav./sum(Expwav);
% Frequency at which the neural data should be sampled
Fs = 1/(TR.*10^-3);
% Loop through the stimuli and fill in the matrix
for stim=1:length(Duration)
    % Time slots for the neural response
    TimeBinsY = -Delay : (Delay + Duration(stim));
    SpikePattern = zeros(1,length(TimeBinsY)-1);
    for isp = 1:length(SAT{stim})
        SpikeInd = round(SAT{stim}(isp));
        if (SpikeInd>=-Delay) && (SpikeInd<(Delay + Duration(stim)))
            SpikePattern(SpikeInd + Delay +1) = SpikePattern(SpikeInd + Delay +1) +1;
        end
    end
    if Fs~=1000
        YPerStim{stim} = resample(conv(SpikePattern, Expwav,'same'), Fs, 1000);
    else
        YPerStim{stim} = conv(SpikePattern, Expwav,'same');
    end
%     % change zero values for the smallest value under matlab.
%     if sum(YPerStim{stim}==0)
%         MinData = min(YPerStim{stim}(YPerStim{stim} ~=0));
%         if ~isempty(MinData)
%             YPerStim{stim}(YPerStim{stim}==0)=min(MinData,realmin('double'));
%         else
%             YPerStim{stim}(YPerStim{stim}==0)=realmin('double');
%         end
%     end
    
    % Make sure that the output sum(Y) = input sum(Y)
    if (~sum(YPerStim{stim})==0) && ~(sum(SpikePattern)==0)
        YPerStim{stim} = YPerStim{stim}./sum(YPerStim{stim}).*sum(SpikePattern);
    end
    
    
    %         % Time slots for the neural response
    %         TimeBinsY = -(Delay) : TR: (Delay + Duration(stim));
    %         YPerStim{stim} = nan(1,length(TimeBinsY)-1);
    %         for tt=1:(length(TimeBinsY)-1)
    %             % Find the number of spikes
    %             YPerStim{stim}(tt) = sum( (SAT{stim}>=TimeBinsY(tt)) .* (SAT{stim}<TimeBinsY(tt+1)));
    %         end
end

end

%     % Calculate the acoustic feature
%     % Temporal Enveloppe
%     % bandpass filter the ambient mic recording
%     Filt_RawVoc = filtfilt(sos_raw_band,1,D2.SpikeTimesVoc.Logger16.VocWave{stim});
%     Amp_env_Mic{stim} = cell(length(Fhigh_powers),1);
%     Power_env_Mic{stim} = cell(length(Fhigh_powers),1);
%     for ff=1:length(Fhigh_powers)
%         [Amp_env_Mic{stim}{ff}, Power_env_Mic{stim}{ff}] = running_rms(Filt_RawVoc, D1.FS, Fhigh_powers(ff), Fs_env);
%     end

function plotxyfeatures(BioSound,YPerStim,XPerStim,TR,Win,Delay,Duration,F_high,FeatureName)
% This function plots for each stimulus the spectrogram with the
% corresponding acoustic feature and spike rate
if nargin<6
    FeatureName = 'Acoustic Feature';
end
DBNOISE =60;
f_low = 0;
close all
for stim =1:length(BioSound)
    figure(1)
    clf
    ColorCode = get(groot,'DefaultAxesColorOrder');
    subplot(2,1,1)
    title(sprintf('vocalization %d/%d', stim,length(BioSound)))
    yyaxis left
    logB = BioSound{stim}.spectro;
    maxB = max(max(logB));
    minB = maxB-DBNOISE;
    imagesc(double(BioSound{stim}.to)*1000,double(BioSound{stim}.fo),logB);          % to is in seconds
    axis xy;
    caxis('manual');
    caxis([minB maxB]);
    cmap = spec_cmap();
    colormap(cmap);
    %         colorbar()
    v_axis = axis;
    v_axis(3)=f_low;
    v_axis(4)=F_high;
    axis(v_axis);
    xlabel('time (ms)'), ylabel('Frequency');
    
    
    hold on
    yyaxis right
    StimFeature = [XPerStim{stim}(:,1)' XPerStim{stim}(end,2:end)];
    MaxT = floor(max(double(BioSound{stim}.to)*1000));
    XStimFeature = [-Win:-1 double(BioSound{stim}.to)*1000 MaxT+(1:(Win))];
    plot(XStimFeature(1:length(StimFeature)), StimFeature, 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    %     IndMax = floor(max(double(BioSound{stim}.to)*1000)/TR);
    %     plot(TR/2 + (0:TR:((IndMax-1)*TR)), StimFeature(1:IndMax), 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    ylabel(sprintf('%s', FeatureName))
    if strcmp(FeatureName,'Spectral Mean')
        ylim(v_axis(3:4))
    end
    xlim([-Win Win+Duration(stim)])
    hold off
    
    subplot(2,1,2)
    yyaxis left
    bar(-(Win-Delay) : (Delay + Duration(stim)-1),YPerStim{stim})
    %     bar(TR/2+(-(Delay) : TR : (Delay + Duration(stim)-TR)),YPerStim{stim})
    ylabel(sprintf('Spikes/ms (%d ms Gauss win)', TR))
    ylim([0 1])
    xlabel('Time (ms)')
    hold on
    line([-(Win-Delay) Duration(stim)+Delay], [0.975 0.975], 'Color',[0 0.4470 0.7410], 'LineWidth',12)
    hold on
    text(0,0.975,'Neural response', 'Color', [1 1 1])
    hold on
    yyaxis right
    plot(XStimFeature(1:length(StimFeature)), StimFeature, 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    %     plot(TR/2 + (0:TR:((IndMax-1)*TR)), StimFeature(1:IndMax), 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    ylabel(sprintf('%s', FeatureName))
    xlim([-Win Win+Duration(stim)])
    hold off
    pause(1)
end
end


function [MSEVal, YPredict] = find_optimalTR(X,Y,Xval,Yval)
    % First find the optimal hyper parameter Lambda using by default a 5 fold
    % cross-validation optimization on MSE with a ridge regularization
    [~,FitInfo,~] = fitrlinear(X,Y,...
    'OptimizeHyperparameters','Lambda','HyperparameterOptimizationOptions',struct('MaxObjectiveEvaluations',20,'UseParallel',true, 'ShowPlots', true),'Regularization', 'ridge', 'Learner','leastsquares');
    % Second, determine the model Beta values using the optimize Lambda and all
    % the dataset
    [Mdl,~] = fitrlinear(X,Y,...
    'Lambda',FitInfo.Lambda,'Regularization', 'ridge', 'Learner','leastsquares');
    % third, find the fit performance of that model on the Validation Dataset
    % Yval
    MSEVal = loss(Mdl, Xval, Yval);
    YPredict = predict(Mdl,Xval);
    figure(10)
    scatter(Yval,YPredict, 'filled')
    xlabel('Observed log Spkie Rate')
    ylabel('Predicted log Spike RAte')
    YL = ylim;
    XL = xlim;
    ylim([min(YL(1),XL(1)) max(YL(2),XL(2))])
    xlim([min(YL(1),XL(1)) max(YL(2),XL(2))])
    F1 = figure(1);
    close(F1)
    F2 = figure(2);
    close(F2)
end


