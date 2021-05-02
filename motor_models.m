addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'));
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'));
addpath(genpath('/Users/elie/Documents/CODE/LMC'));
addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'));
DatFig=0; %Set to 1 to see input data figures for each cell
OutFig = 1;%Set to 1 to see output data figures for each cell

%% Listing datacells
%Filename = '59834_20190611_SSS_1-97.mat';
% Filename = '59834_20190610_SSS_1-130.mat';
%Path = '/Users/elie/Documents/LMCResults/';
Path = '/Volumes/JulieE8T/LMCResults/';
% Path = '/Users/elie/Documents/ManipBats/LMC/ResultsFiles/';
% Path = '/Users/elie/Google Drive/BatmanData/';

AllFiles = dir(fullfile(Path,'59834*.mat'));
Files2run = zeros(length(AllFiles),1);
for ff=1:length(AllFiles)
    fprintf(1, 'Cell %d/%d\n', ff, length(AllFiles))
    Cell = load(fullfile(AllFiles(ff).folder, AllFiles(ff).name), 'What');
    if (length(strfind(AllFiles(ff).name, '_'))==3) && isfield(Cell, 'What') % Cells without the what field are non-compiled cells probably because they are just noise
        Files2run(ff) = 1;
    end
    clear Cell
end
CellsPath = AllFiles(logical(Files2run));


%% Check cell stability between 2 sessions
NCells = length(CellsPath);
KDE_OpvsFr_p = nan(NCells,1);
KDE_OpvsFr_KStestp = nan(NCells,1);
for cc=1:NCells
    
    %% load data
    CellPath = fullfile(CellsPath(cc).folder,CellsPath(cc).name);
    fprintf(1, '*** Cell %s %d/%d ***\n', CellPath, cc, NCells)
    Cell = load(CellPath,'QualitySSU', 'FreeBehavSession', 'OperantSession');
    % decide if cell is stable enough to compare free and operant sessions
    IndicesOp=((Cell.QualitySSU.TimePoints(2:end)-60/2)>0) .* ((Cell.QualitySSU.TimePoints(2:end)-60/2)<diff(Cell.OperantSession));
    IndicesFr=((Cell.QualitySSU.TimePoints(2:end)-60/2)>(Cell.FreeBehavSession(1)- Cell.OperantSession(1))) .* ((Cell.QualitySSU.TimePoints(2:end)-60/2)<(Cell.FreeBehavSession(2)- Cell.OperantSession(1)));
    [H,KDE_OpvsFr_p(cc),CI,STATS] = ttest2(Cell.QualitySSU.KDE(logical(IndicesOp)), Cell.QualitySSU.KDE(logical(IndicesFr)));
    if sum(IndicesFr)&& sum(IndicesOp)
        [h,KDE_OpvsFr_KStestp(cc),ks2stat] = kstest2(Cell.QualitySSU.KDE(logical(IndicesOp)), Cell.QualitySSU.KDE(logical(IndicesFr)));
    end
%     figure(15);clf;histogram(Cell.QualitySSU.KDE(logical(IndicesOp)), 'Normalization', 'probability'); hold on; histogram(Cell.QualitySSU.KDE(logical(IndicesFr)), 'Normalization', 'probability')
%     xlabel('Spike Rate in Hz (bin=1min)')
%     ylabel('probability')
%     legend({'operant' 'Free'})
%     title(sprintf('Spike rate distribution p=%.2f',KDE_OpvsFr_p(cc)))
%     pause(1)
end  
save(fullfile(Path, 'KDE_SessionStability.mat'), 'CellsPath', 'KDE_OpvsFr_KStestp','KDE_OpvsFr_p')

%% Running through cells to find the optimal time resolution of the neural response for acoustic feature predicion from the neural response
PlotCoherenceFig = 0; % To plot the result of coherence calculation for each cell
%Trill=0; % Set to 1 to only do calculations on Trill calls, set to 0 to do calculations on allcalls
% here we calculate the coherency between the neural response and the
% acoustic features
StimXYDataPlot=0;
FeatureName = 'amp';% {'amp' 'SpectralMean' 'sal'};
SessionChoices = {'Operant', 'Free'};
TR=2; % 2ms is chosen as the Time resolution for the neural data
Fs = 1/(TR*10^-3); % the data is then sampled at the optimal frequency given the neural time resolution choosen
% find the closest power of 2 for the number of FFT window points that
% correspond to the Nyquist limit
Nyquist = Fs * 0.5;
nFFT = 2^ceil(log2(Nyquist));
NCells = length(CellsPath);
%Delay = nFFT/(2*Fs)*10^3;
Delay=200; % The segment of data taken into account is -Delay ms before the vocalization onset and +200ms after the vocalization offset
NBoot = 500; % number of voc ID permutation bootstraps for the significance of Info on coherence for each cell
Trill=0; % 0: all calls; 1: Trills only; 2: non Trills only(still have to write the code for that!!)
BootstrapType=0;% 0: bootstrap calculation based on Voc ID shuffling, keeping onsets; any values other than 0, shuffling in time within +/-BootstrapType ms keeping voc ID intact
NumVoc = 0;% 0: keep all vocalizations; any value: randomly choose NumVoc vocalizations to calculate Coherence
% Lags = -Delay:Delay;
% Freqs = (0:ceil(length(Lags)/2)).* (2*Nyquist/length(Lags)); % Lags is a uneven number so F(i) = i*2*Nyquist/length(Lags)
AuditoryCoherenceFreeAll.CoherencyT_DelayAtzero = nan(NCells,1);
AuditoryCoherenceFreeAll.CoherencyT_WidthAtMaxPeak = nan(NCells,1);
AuditoryCoherenceFreeAll.MaxCoherence = nan(NCells,1);
AuditoryCoherenceFreeAll.MaxCoherenceF = nan(NCells,1);
AuditoryCoherenceFreeAll.CoherencePeaks = cell(NCells,1);
AuditoryCoherenceFreeAll.CoherencePeaksF = cell(NCells,1);
AuditoryCoherenceFreeAll.FirstNonSigCoherenceFreq = nan(NCells,1);
AuditoryCoherenceFreeAll.SecondCoherenceFreqCutOff = nan(NCells,1);
AuditoryCoherenceFreeAll.Info = nan(NCells,1);
AuditoryCoherenceFreeAll.Info_low = nan(NCells,1);
AuditoryCoherenceFreeAll.Info_up = nan(NCells,1);
AuditoryCoherenceFreeAll.Info_p = nan(NCells,1);
AuditoryCoherenceFreeAll.CellsPath = CellsPath;
AuditoryCoherenceFreeAll.TR = TR;
AuditoryCoherenceFreeAll.Delay = Delay;

AuditoryCoherenceOperantAll = AuditoryCoherenceFreeAll;
% MotorCoherenceFreeAll = AuditoryCoherenceFreeAll;
% MotorCoherenceOperantAll = AuditoryCoherenceFreeAll;

for cc=1:NCells
    
    CellTimer = tic();
    %% load data
    CellPath = fullfile(CellsPath(cc).folder,CellsPath(cc).name);
    fprintf(1, '*** Cell %s %d/%d %s ***\n', CellPath, cc, NCells, FeatureName)
    Cell = load(CellPath, 'What', 'ExpType', 'Who', 'BioSound', 'Duration','SpikesArrivalTimes_Behav','QualitySSU', 'VocOverlap', 'AudioQuality', 'DelayBefore','DelayAfter');
    
    if ~Trill
        
        AuditoryCoherenceFree = coherence4LMC(Cell,0,'Free', 'amp', Trill,NumVoc,BootstrapType);
        
%         MotorCoherenceFree = coherence4LMC(Cell, 1,'Free', 'amp', Trill,NumVoc,BootstrapType);
        
        AuditoryCoherenceOperant = coherence4LMC(Cell, 0,'Operant', 'amp', Trill,NumVoc,BootstrapType);
        
%         MotorCoherenceOperant = coherence4LMC(Cell, 1,'Operant', 'amp', Trill,NumVoc,BootstrapType);
        if length(BootstrapType)==1 && (BootstrapType==0) % only save to cell file if we are on the VocID shuffling for now
            save(CellPath, 'AuditoryCoherenceFree','AuditoryCoherenceOperant', '-append')
%             save(CellPath, 'AuditoryCoherenceFree','MotorCoherenceFree','AuditoryCoherenceOperant', 'MotorCoherenceOperant','-append')
        end
        if ~isfield(AuditoryCoherenceFree, 'Error')
            AuditoryCoherenceFreeAll.CoherencyT_DelayAtzero(cc) = AuditoryCoherenceFree.CoherencyT_DelayAtzero;
            AuditoryCoherenceFreeAll.CoherencyT_WidthAtMaxPeak(cc) = AuditoryCoherenceFree.CoherencyT_WidthAtMaxPeak;
            AuditoryCoherenceFreeAll.MaxCoherence(cc) = AuditoryCoherenceFree.MaxCoherence1;
            AuditoryCoherenceFreeAll.MaxCoherenceF(cc) = AuditoryCoherenceFree.MaxCoherence2;
            AuditoryCoherenceFreeAll.CoherencePeaks{cc} = AuditoryCoherenceFree.CoherencePeaks;
            AuditoryCoherenceFreeAll.CoherencePeaksF{cc} = AuditoryCoherenceFree.CoherencePeaksF;
            AuditoryCoherenceFreeAll.FirstNonSigCoherenceFreq(cc) = AuditoryCoherenceFree.FirstNonSigCoherenceFreq;
            AuditoryCoherenceFreeAll.SecondCoherenceFreqCutOff(cc) = AuditoryCoherenceFree.SecondCoherenceFreqCutOff;
            AuditoryCoherenceFreeAll.Info(cc) = AuditoryCoherenceFree.Info;
            AuditoryCoherenceFreeAll.Info_low(cc) = AuditoryCoherenceFree.Info_low;
            AuditoryCoherenceFreeAll.Info_up(cc) = AuditoryCoherenceFree.Info_up;
            AuditoryCoherenceFreeAll.Info_p(cc) = AuditoryCoherenceFree.Info_p;
        end
%         if ~isfield(MotorCoherenceFree, 'Error')
%             MotorCoherenceFreeAll.CoherencyT_DelayAtzero(cc) = MotorCoherenceFree.CoherencyT_DelayAtzero;
%             MotorCoherenceFreeAll.CoherencyT_WidthAtMaxPeak(cc) = MotorCoherenceFree.CoherencyT_WidthAtMaxPeak;
%             MotorCoherenceFreeAll.MaxCoherence(cc) = MotorCoherenceFree.MaxCoherence1;
%             MotorCoherenceFreeAll.MaxCoherenceF(cc) = MotorCoherenceFree.MaxCoherence2;
%             MotorCoherenceFreeAll.CoherencePeaks{cc} = MotorCoherenceFree.CoherencePeaks;
%             MotorCoherenceFreeAll.CoherencePeaksF{cc} = MotorCoherenceFree.CoherencePeaksF;
%             MotorCoherenceFreeAll.FirstNonSigCoherenceFreq(cc) = MotorCoherenceFree.FirstNonSigCoherenceFreq;
%             MotorCoherenceFreeAll.SecondCoherenceFreqCutOff(cc) = MotorCoherenceFree.SecondCoherenceFreqCutOff;
%             MotorCoherenceFreeAll.Info(cc) = MotorCoherenceFree.Info;
%             MotorCoherenceFreeAll.Info_low(cc) = MotorCoherenceFree.Info_low;
%             MotorCoherenceFreeAll.Info_up(cc) = MotorCoherenceFree.Info_up;
%             MotorCoherenceFreeAll.Info_p(cc) = MotorCoherenceFree.Info_p;
%         end

        if ~isfield(AuditoryCoherenceOperant, 'Error')
            AuditoryCoherenceOperantAll.CoherencyT_DelayAtzero(cc) = AuditoryCoherenceOperant.CoherencyT_DelayAtzero;
            AuditoryCoherenceOperantAll.CoherencyT_WidthAtMaxPeak(cc) = AuditoryCoherenceOperant.CoherencyT_WidthAtMaxPeak;
            AuditoryCoherenceOperantAll.MaxCoherence(cc) = AuditoryCoherenceOperant.MaxCoherence1;
            AuditoryCoherenceOperantAll.MaxCoherenceF(cc) = AuditoryCoherenceOperant.MaxCoherence2;
            AuditoryCoherenceOperantAll.CoherencePeaks{cc} = AuditoryCoherenceOperant.CoherencePeaks;
            AuditoryCoherenceOperantAll.CoherencePeaksF{cc} = AuditoryCoherenceOperant.CoherencePeaksF;
            AuditoryCoherenceOperantAll.FirstNonSigCoherenceFreq(cc) = AuditoryCoherenceOperant.FirstNonSigCoherenceFreq;
            AuditoryCoherenceOperantAll.SecondCoherenceFreqCutOff(cc) = AuditoryCoherenceOperant.SecondCoherenceFreqCutOff;
            AuditoryCoherenceOperantAll.Info(cc) = AuditoryCoherenceOperant.Info;
            AuditoryCoherenceOperantAll.Info_low(cc) = AuditoryCoherenceOperant.Info_low;
            AuditoryCoherenceOperantAll.Info_up(cc) = AuditoryCoherenceOperant.Info_up;
            AuditoryCoherenceOperantAll.Info_p(cc) = AuditoryCoherenceOperant.Info_p;
        end

%         if ~isfield(MotorCoherenceOperant, 'Error')
%             MotorCoherenceOperantAll.CoherencyT_DelayAtzero(cc) = MotorCoherenceOperant.CoherencyT_DelayAtzero;
%             MotorCoherenceOperantAll.CoherencyT_WidthAtMaxPeak(cc) = MotorCoherenceOperant.CoherencyT_WidthAtMaxPeak;
%             MotorCoherenceOperantAll.MaxCoherence(cc) = MotorCoherenceOperant.MaxCoherence1;
%             MotorCoherenceOperantAll.MaxCoherenceF(cc) = MotorCoherenceOperant.MaxCoherence2;
%             MotorCoherenceOperantAll.CoherencePeaks{cc} = MotorCoherenceOperant.CoherencePeaks;
%             MotorCoherenceOperantAll.CoherencePeaksF{cc} = MotorCoherenceOperant.CoherencePeaksF;
%             MotorCoherenceOperantAll.FirstNonSigCoherenceFreq(cc) = MotorCoherenceOperant.FirstNonSigCoherenceFreq;
%             MotorCoherenceOperantAll.SecondCoherenceFreqCutOff(cc) = MotorCoherenceOperant.SecondCoherenceFreqCutOff;
%             MotorCoherenceOperantAll.Info(cc) = MotorCoherenceOperant.Info;
%             MotorCoherenceOperantAll.Info_low(cc) = MotorCoherenceOperant.Info_low;
%             MotorCoherenceOperantAll.Info_up(cc) = MotorCoherenceOperant.Info_up;
%             MotorCoherenceOperantAll.Info_p(cc) = MotorCoherenceOperant.Info_p;
%         end
    elseif Trill==1 % Only calculating motor coherence for Trills because only Cooper is producing trills and never listening to them
        MotorCoherenceFreeTrill = coherence4LMC(Cell, 1,'Free', 'amp', Trill,NumVoc, BootstrapType);

        MotorCoherenceOperantTrill = coherence4LMC(Cell, 1,'Operant', 'amp', Trill,NumVoc, BootstrapType);
        if length(BootstrapType)==1 && (BootstrapType(1)==0) % only save to cell file if we are on the VocID shuffling for now
            save(CellPath, 'MotorCoherenceFreeTrill', 'MotorCoherenceOperantTrill','-append')
        end

        if ~isfield(MotorCoherenceFreeTrill, 'Error')
            MotorCoherenceFreeAll.CoherencyT_DelayAtzero(cc) = MotorCoherenceFreeTrill.CoherencyT_DelayAtzero;
            MotorCoherenceFreeAll.CoherencyT_WidthAtMaxPeak(cc) = MotorCoherenceFreeTrill.CoherencyT_WidthAtMaxPeak;
            MotorCoherenceFreeAll.MaxCoherence(cc) = MotorCoherenceFreeTrill.MaxCoherence1;
            MotorCoherenceFreeAll.MaxCoherenceF(cc) = MotorCoherenceFreeTrill.MaxCoherence2;
            MotorCoherenceFreeAll.CoherencePeaks{cc} = MotorCoherenceFreeTrill.CoherencePeaks;
            MotorCoherenceFreeAll.CoherencePeaksF{cc} = MotorCoherenceFreeTrill.CoherencePeaksF;
            MotorCoherenceFreeAll.FirstNonSigCoherenceFreq(cc) = MotorCoherenceFreeTrill.FirstNonSigCoherenceFreq;
            MotorCoherenceFreeAll.SecondCoherenceFreqCutOff(cc) = MotorCoherenceFreeTrill.SecondCoherenceFreqCutOff;
            MotorCoherenceFreeAll.Info(cc) = MotorCoherenceFreeTrill.Info;
            MotorCoherenceFreeAll.Info_low(cc) = MotorCoherenceFreeTrill.Info_low;
            MotorCoherenceFreeAll.Info_up(cc) = MotorCoherenceFreeTrill.Info_up;
            MotorCoherenceFreeAll.Info_p(cc) = MotorCoherenceFreeTrill.Info_p;
        end


        if ~isfield(MotorCoherenceOperantTrill, 'Error')
            MotorCoherenceOperantAll.CoherencyT_DelayAtzero(cc) = MotorCoherenceOperantTrill.CoherencyT_DelayAtzero;
            MotorCoherenceOperantAll.CoherencyT_WidthAtMaxPeak(cc) = MotorCoherenceOperantTrill.CoherencyT_WidthAtMaxPeak;
            MotorCoherenceOperantAll.MaxCoherence(cc) = MotorCoherenceOperantTrill.MaxCoherence1;
            MotorCoherenceOperantAll.MaxCoherenceF(cc) = MotorCoherenceOperantTrill.MaxCoherence2;
            MotorCoherenceOperantAll.CoherencePeaks{cc} = MotorCoherenceOperantTrill.CoherencePeaks;
            MotorCoherenceOperantAll.CoherencePeaksF{cc} = MotorCoherenceOperantTrill.CoherencePeaksF;
            MotorCoherenceOperantAll.FirstNonSigCoherenceFreq(cc) = MotorCoherenceOperantTrill.FirstNonSigCoherenceFreq;
            MotorCoherenceOperantAll.SecondCoherenceFreqCutOff(cc) = MotorCoherenceOperantTrill.SecondCoherenceFreqCutOff;
            MotorCoherenceOperantAll.Info(cc) = MotorCoherenceOperantTrill.Info;
            MotorCoherenceOperantAll.Info_low(cc) = MotorCoherenceOperantTrill.Info_low;
            MotorCoherenceOperantAll.Info_up(cc) = MotorCoherenceOperantTrill.Info_up;
            MotorCoherenceOperantAll.Info_p(cc) = MotorCoherenceOperantTrill.Info_p;
        end
    end
        
        
    fprintf(1, '******************** Done with Cell %s %d/%d in %ds ******************\n', CellPath,cc,NCells, toc(CellTimer));
    %             fprintf(1, 'Done with Cell %d/%d\n', cc,NCells)
    clear Cell
end

% Order cells by
% decreasing values of info
if ~Trill && length(BootstrapType)==1 && (BootstrapType(1)==0)
    [~,AuditoryCoherenceFreeAll.GoodInfo] = sort(AuditoryCoherenceFreeAll.Info,'descend');
    save(fullfile(Path,sprintf('AuditoryCoherence_%s_%s.mat', FeatureName, 'Free')),'-struct','AuditoryCoherenceFreeAll');
    
%     [~,MotorCoherenceFreeAll.GoodInfo] = sort(MotorCoherenceFreeAll.Info,'descend');
%     save(fullfile(Path,sprintf('MotorCoherence_%s_%s.mat', FeatureName, 'Free')),'-struct','MotorCoherenceFreeAll');
%     
    [~,AuditoryCoherenceOperantAll.GoodInfo] = sort(AuditoryCoherenceOperantAll.Info,'descend');
    save(fullfile(Path,sprintf('AuditoryCoherence_%s_%s.mat', FeatureName, 'Operant')),'-struct','AuditoryCoherenceOperantAll');
    
%     [~,MotorCoherenceOperantAll.GoodInfo] = sort(MotorCoherenceOperantAll.Info,'descend');
%     save(fullfile(Path,sprintf('MotorCoherence_%s_%s.mat', FeatureName, 'Operant')),'-struct','MotorCoherenceOperantAll');
elseif Trill && length(BootstrapType)==1 && (BootstrapType(1)==0)
    [~,MotorCoherenceFreeAll.GoodInfo] = sort(MotorCoherenceFreeAll.Info,'descend');
    save(fullfile(Path,sprintf('MotorCoherence_%s_%s_Trill.mat', FeatureName, 'Free')),'-struct','MotorCoherenceFreeAll');
    
    [~,MotorCoherenceOperantAll.GoodInfo] = sort(MotorCoherenceOperantAll.Info,'descend');
    save(fullfile(Path,sprintf('MotorCoherence_%s_%s_Trill.mat', FeatureName, 'Operant')),'-struct','MotorCoherenceOperantAll');
elseif Trill && (~length(BootstrapType)==1 || (BootstrapType(1)~=0))
    [~,MotorCoherenceFreeAll.GoodInfo] = sort(MotorCoherenceFreeAll.Info,'descend');
    if length(BootstrapType)==1
        save(fullfile(Path,sprintf('MotorCoherence_%s_%s_Trill_TimeBootstrap%d.mat', FeatureName, 'Free',BootstrapType(1))),'-struct','MotorCoherenceFreeAll');
    elseif length(BootstrapType)==2
        save(fullfile(Path,sprintf('MotorCoherence_%s_%s_Trill_TimeBootstrap%d-%d.mat', FeatureName, 'Free',BootstrapType(1),BootstrapType(2))),'-struct','MotorCoherenceFreeAll');
    end
    [~,MotorCoherenceOperantAll.GoodInfo] = sort(MotorCoherenceOperantAll.Info,'descend');
    if length(BootstrapType)==1
        save(fullfile(Path,sprintf('MotorCoherence_%s_%s_Trill_TimeBootstrap%d.mat', FeatureName, 'Operant',BootstrapType)),'-struct','MotorCoherenceOperantAll');
    elseif length(BootstrapType)==2
        save(fullfile(Path,sprintf('MotorCoherence_%s_%s_Trill_TimeBootstrap%d-%d.mat', FeatureName, 'Operant',BootstrapType(1),BootstrapType(2))),'-struct','MotorCoherenceOperantAll');
    end
elseif ~Trill && (length(BootstrapType)==1 && (BootstrapType(1)~=0))
    [~,AuditoryCoherenceFreeAll.GoodInfo] = sort(AuditoryCoherenceFreeAll.Info,'descend');
    save(fullfile(Path,sprintf('AuditoryCoherence_%s_%s_TimeBootstrap%d.mat', FeatureName, 'Free',BootstrapType)),'-struct','AuditoryCoherenceFreeAll');
    
    [~,MotorCoherenceFreeAll.GoodInfo] = sort(MotorCoherenceFreeAll.Info,'descend');
    save(fullfile(Path,sprintf('MotorCoherence_%s_%s_TimeBootstrap%d.mat', FeatureName, 'Free',BootstrapType)),'-struct','MotorCoherenceFreeAll');
    
    [~,AuditoryCoherenceOperantAll.GoodInfo] = sort(AuditoryCoherenceOperantAll.Info,'descend');
    save(fullfile(Path,sprintf('AuditoryCoherence_%s_%s_TimeBootstrap%d.mat', FeatureName, 'Operant',BootstrapType)),'-struct','AuditoryCoherenceOperantAll');
    
    [~,MotorCoherenceOperantAll.GoodInfo] = sort(MotorCoherenceOperantAll.Info,'descend');
    save(fullfile(Path,sprintf('MotorCoherence_%s_%s_TimeBootstrap%d.mat', FeatureName, 'Operant',BootstrapType)),'-struct','MotorCoherenceOperantAll');
else
    warning('Code not written for these conditions!')
    keyboard
end

%% Plots results of coherence calculations for the population
FeatureName = {'amp' 'SpectralMean' 'sal'};
SessionChoice = 'Free'; % Operant or Free
MS = 'Auditory'; % Motor or Auditory
Trill=0;
% for fn = 1:length(FeatureName)
    fn=1;
    if ~Trill
        load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, SessionChoice)))
    elseif Trill==1
        load(fullfile(Path,sprintf('%sCoherence_%s_%s_Trill.mat', MS, FeatureName{fn}, SessionChoice)))
    end
    if strcmp(FeatureName{fn}, 'amp')
        FeatureName2 = 'Amplitude';
    elseif strcmp(FeatureName{fn}, 'sal')
        FeatureName2 = 'Pitch Saliency';
    elseif strcmp(FeatureName{fn}, 'SpectralMean')
        FeatureName2 = 'Spectral Mean';
    end
    figure(1)
    clf
    subplot(1,3,1)
    histogram(MaxCoherence(:,1),'BinWidth',0.005,'FaceColor','k')
    xlabel(sprintf('Max coherence with sound %s', FeatureName2))
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
    XLIM = get(gca, 'XLim');
    set(gca, 'XLim', [XLIM(1) 0.55])
%     keyboard
    
    subplot(1,3,2)
    Colors = [(Info_p<0.01) zeros(length(Info),1) zeros(length(Info),1)];
    scatter(Info,MaxCoherence(:,1),20, Colors, 'filled')
%     plot(Info,MaxCoherence(:,1), 'o','Color','k','MarkerSize',6,'MarkerFaceColor','k')
    xlabel(sprintf('Information on Coherence with sound %s (bits)', FeatureName2))
    ylabel(sprintf('Max coherence with sound %s',FeatureName2))
    YLIM = get(gca, 'YLim');
    set(gca, 'YLim', [YLIM(1) 0.6])
    XLIM = get(gca, 'XLim');
    set(gca, 'XLim', [XLIM(1) 5])
%     keyboard
    
    subplot(1,3,3)
    histogram(Info,'BinWidth',0.05,'FaceColor','k')
    xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
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
    XLIM = get(gca, 'XLim');
    set(gca, 'XLim', [XLIM(1) 5])
%     keyboard
if ~Trill
    suplabel(sprintf('%s Coherence %s %s session', MS, FeatureName2, SessionChoice), 't');
elseif Trill==1
    suplabel(sprintf('%s Coherence %s Trills only %s session', MS, FeatureName2, SessionChoice), 't');
end
    
    figure(2)
    clf
    subplot(2,5,1)
    Colors = [(Info_p<0.01) zeros(length(Info),1) zeros(length(Info),1)];
    scatter(Info,CoherencyT_DelayAtzero,20, Colors, 'filled')
    xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
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
    scatter(Info,CoherencyT_WidthAtMaxPeak,20, Colors, 'filled')
    xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
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
    scatter(Info,FirstNonSigCoherenceFreq,20, Colors, 'filled')
    xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
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
    scatter(Info,SecondCoherenceFreqCutOff,20, Colors, 'filled')
    xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
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
    scatter(Info,MaxCoherenceF,20, Colors, 'filled')
    xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
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
    histogram(CoherencyT_DelayAtzero(Info_p>=0.01),'BinWidth',TR/2,'FaceColor','k')
    hold on
    histogram(CoherencyT_DelayAtzero(Info_p<0.01),'BinWidth',TR/2,'FaceColor','r')
    xlabel('Phase of coherency (ms)')
    ylabel('Number of cells')
    hold on
    v=vline(0, 'r:');
    v.LineWidth = 2;
    hold on
    text(-57, 27, sprintf('All mean = %.2f +/- %.2f', nanmean(CoherencyT_DelayAtzero), nanstd(CoherencyT_DelayAtzero)/(sum(~isnan(CoherencyT_DelayAtzero)))^0.5))
    hold on
    text(-57, 29, sprintf('Sig mean = %.2f +/- %.2f', nanmean(CoherencyT_DelayAtzero(Info_p<0.01)), nanstd(CoherencyT_DelayAtzero(Info_p<0.01))/(sum(~isnan(CoherencyT_DelayAtzero(Info_p<0.01))))^0.5))
    
    subplot(2,5,7)
    histogram(CoherencyT_WidthAtMaxPeak(Info_p>=0.01),'BinWidth',TR/2,'FaceColor','k')
    hold on
    histogram(CoherencyT_WidthAtMaxPeak(Info_p<0.01),'BinWidth',TR/2,'FaceColor','r')
    xlabel('Time resolution of Coherency (ms)')
    ylabel('Number of cells')
    hold on
    text(35, 32, sprintf('All mean = %.2f +/- %.2f', nanmean(CoherencyT_WidthAtMaxPeak), nanstd(CoherencyT_WidthAtMaxPeak)/(sum(~isnan(CoherencyT_WidthAtMaxPeak)))^0.5))
    hold on
    text(35, 34, sprintf('Sig mean = %.2f +/- %.2f', nanmean(CoherencyT_WidthAtMaxPeak(Info_p<0.01)), nanstd(CoherencyT_WidthAtMaxPeak(Info_p<0.01))/(sum(~isnan(CoherencyT_WidthAtMaxPeak(Info_p<0.01))))^0.5))
    
    
    
    subplot(2,5,8)
    histogram(FirstNonSigCoherenceFreq(Info_p>=0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
    hold on
    histogram(FirstNonSigCoherenceFreq(Info_p<0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','r')
    xlabel('Max significant Frequency (Hz)')
    ylabel('Number of Cells')
    
    subplot(2,5,9)
    histogram(SecondCoherenceFreqCutOff(Info_p>=0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
    hold on
    histogram(SecondCoherenceFreqCutOff(Info_p<0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','r')
    xlabel('Max 2nd peaks significant Frequency (Hz)')
    ylabel('Number of Cells')
    
    subplot(2,5,10)
    histogram(MaxCoherenceF(Info_p>=0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
    hold on
    histogram(MaxCoherenceF(Info_p<0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','r')
    % Some points with very low values of info have high values of frequency of
    % Max coherence, keep the plot focused on the majority of points
    xlim([0 50])
    xlabel('Frequency of Max Coherence (Hz)')
    ylabel('Number of Cells')
    
    if ~Trill
        suplabel(sprintf('%s Coherence %s session %s N=%d cells %d with permutation test p<0.01', MS, SessionChoice, FeatureName2, sum(~isnan(Info)), sum(Info_p<0.01)), 't');
    elseif Trill==1
        suplabel(sprintf('%s Coherence %s Trills only session %s N=%d cells %d with permutation test p<0.01', MS, SessionChoice, FeatureName2, sum(~isnan(Info)), sum(Info_p<0.01)), 't');
    end
    
    fprintf(1,'Cell with highest Info Value: %s', CellsPath(Info == max(Info)).name)
    % Amplitude: 59834_20190614_SSS_1-100.mat (cc=91/488)
    % SpectralMean 59834_20190610_SSS_1-130.mat (cc=43/488)
    % Saliency: 59834_20190708_SSM_1-228.mat (cc=454/488)
    
    % Plot the values of the secondary peaks found for some cells
    figure(4)
    scatter([CoherencePeaksF{:}], [CoherencePeaks{:}],10, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor','k')
    xlabel('Frequency of secondary peaks in Coherence')
    ylabel('Values of Coherence')
    
    
    
    
    %     % plot again the histograms for high Info subset of cells
    %     figure(5)
    %     clf
    %     subplot(1,5,1)
    %     histogram(CoherencyT_DelayAtzero(GoodInfo),'BinWidth',TR/2,'FaceColor','k')
    %     xlabel('Phase of coherency (ms)')
    %     ylabel('Number of cells with Info > 1bit')
    %     hold on
    %     v=vline(0, 'r:');
    %     v.LineWidth = 2;
    %
    %     subplot(1,5,2)
    %     histogram(CoherencyT_WidthAtMaxPeak(GoodInfo),'BinWidth',TR/2,'FaceColor','k')
    %     xlabel('Time resolution of Coherency (ms)')
    %     ylabel('Number of cells with Info > 1bit')
    %
    %
    %     subplot(1,5,3)
    %     histogram(FirstNonSigCoherenceFreq(GoodInfo), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
    %     xlabel('Max significant Frequency (Hz)')
    %     ylabel('Number of cells with Info > 1bit')
    %
    %     subplot(1,5,4)
    %     histogram(SecondCoherenceFreqCutOff(GoodInfo), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
    %     xlabel('Max 2nd peaks significant Frequency (Hz)')
    %     ylabel('Number of Cells')
    %
    %     subplot(1,5,5)
    %     histogram(MaxCoherence(find(GoodInfo),2), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
    %     % Some points with very low values of info have high values of frequency of
    %     % Max coherence, keep the plot focused on the majority of points
    %     xlim([0 50])
    %     xlabel('Frequency of Max Coherence (Hz)')
    %     ylabel('Number of cells with Info > 1bit')
    %
    %     suplabel(sprintf('%s',FeatureName2),'t');
%     pause();
% end

%% Comparing values of information, Delay and Time resolution for saliency and amplitude
SessionChoice = 'Operant';

CoherenceSal = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s.mat', 'sal',SessionChoice)));
CoherenceAmp = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s.mat', 'amp',SessionChoice)));
CoherenceSpecMean = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s.mat', 'SpectralMean',SessionChoice)));
figure(6)
clf
subplot(3,3,1)
scatter(CoherenceAmp.Info, CoherenceSal.Info,40,[0 0 0], 'filled')
hold on
plot([0 5], [0 5], 'r:', 'LineWidth',2)
hold off
xlabel('Information on Coherence with sound Amplitude (bits)')
ylabel('Information on Coherence with pitch saliency (bits)')


subplot(3,3,2)
scatter(CoherenceAmp.CoherencyT_DelayAtzero, CoherenceSal.CoherencyT_DelayAtzero,40,[0 0 0], 'filled')
hold on
plot([-60 60], [-60 60], 'r:', 'LineWidth',2)
hold off
xlabel('Phase of Coherency with sound Amplitude (ms)')
ylabel('Phase of Coherency with pitch saliency (ms)')

subplot(3,3,3)
scatter(CoherenceAmp.CoherencyT_WidthAtMaxPeak, CoherenceSal.CoherencyT_WidthAtMaxPeak,40,[0 0 0], 'filled')
hold on
plot([40 350], [50 350], 'r:', 'LineWidth',2)
hold off
xlabel('Time resolution of Coherency with sound Amplitude (ms)')
ylabel('Time resolution of Coherency with pitch saliency (ms)')

subplot(3,3,4)
scatter(CoherenceAmp.Info, CoherenceSpecMean.Info,40,[0 0 0], 'filled')
hold on
plot([0 5], [0 5], 'r:', 'LineWidth',2)
hold off
xlabel('Information on Coherence with sound Amplitude (bits)')
ylabel('Information on Coherence with Spectral Mean (bits)')


subplot(3,3,5)
scatter(CoherenceAmp.CoherencyT_DelayAtzero, CoherenceSpecMean.CoherencyT_DelayAtzero,40,[0 0 0], 'filled')
hold on
plot([-60 60], [-60 60], 'r:', 'LineWidth',2)
hold off
xlabel('Phase of Coherency with sound Amplitude (ms)')
ylabel('Phase of Coherency with Spectral Mean (ms)')

subplot(3,3,6)
scatter(CoherenceAmp.CoherencyT_WidthAtMaxPeak, CoherenceSpecMean.CoherencyT_WidthAtMaxPeak,40,[0 0 0], 'filled')
hold on
plot([40 350], [50 350], 'r:', 'LineWidth',2)
hold off
xlabel('Time resolution of Coherency with sound Amplitude (ms)')
ylabel('Time resolution of Coherency with Spectral Mean (ms)')

subplot(3,3,7)
scatter(CoherenceSal.Info, CoherenceSpecMean.Info,40,[0 0 0], 'filled')
hold on
plot([0 5], [0 5], 'r:', 'LineWidth',2)
hold off
xlabel('Information on Coherence with pitch saliency (bits)')
ylabel('Information on Coherence with Spectral Mean (bits)')


subplot(3,3,8)
scatter(CoherenceSal.CoherencyT_DelayAtzero, CoherenceSpecMean.CoherencyT_DelayAtzero,40,[0 0 0], 'filled')
hold on
plot([-60 60], [-60 60], 'r:', 'LineWidth',2)
hold off
xlabel('Phase of Coherency with pitch saliency (ms)')
ylabel('Phase of Coherency with Spectral Mean (ms)')

subplot(3,3,9)
scatter(CoherenceSal.CoherencyT_WidthAtMaxPeak, CoherenceSpecMean.CoherencyT_WidthAtMaxPeak,40,[0 0 0], 'filled')
hold on
plot([40 350], [50 350], 'r:', 'LineWidth',2)
hold off
xlabel('Time resolution of Coherency with pitch saliency (ms)')
ylabel('Time resolution of Coherency with Spectral Mean (ms)')

suplabel(sprintf('%s: Coherence of Pitch saliency vs Amplitude vs Spectral Mean', SessionChoice), 't')


%% Compare Information on Coherence between Free and Operant session
FeatureName = {'amp' 'SpectralMean' 'sal'};
MS = 'Motor'; % Motor or Auditory
Trill=0;
% for fn = 1:length(FeatureName)
fn=1;
figure(5)
clf
% for fn = 1:length(FeatureName)
    fn=1;
    if ~Trill
        Free = load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, 'Free')));
        Operant = load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, 'Operant')));
        
    elseif Trill==1
        Free = load(fullfile(Path,sprintf('%sCoherence_%s_%s_Trill.mat', MS, FeatureName{fn}, 'Free')));
        Operant = load(fullfile(Path,sprintf('%sCoherence_%s_%s_Trill.mat', MS, FeatureName{fn}, 'Operant')));
        
    end
    Colors = [(Operant.Info_p<0.01) zeros(length(Operant.Info),1) (Free.Info_p<0.01)];
    
    if strcmp(FeatureName{fn}, 'amp')
        FeatureName2 = 'Amplitude';
    elseif strcmp(FeatureName{fn}, 'sal')
        FeatureName2 = 'Pitch Saliency';
    elseif strcmp(FeatureName{fn}, 'SpectralMean')
        FeatureName2 = 'Spectral Mean';
    end
    subplot(1,2,1)
    scatter(Operant.Info,Free.Info,40,Colors, 'filled')
    hold on
    plot([0 7], [0 7], 'r:', 'LineWidth',2)
    hold off
    xlabel('Operant Session')
    ylabel('Free Session')
    ValidCells4Plot = logical((~isnan(Operant.Info .* Free.Info)));
    if ~Trill
        title(sprintf('All Cells (bits; N=%d)', sum(ValidCells4Plot)))
        text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
    elseif Trill==1
        title(sprintf('All Cells (bits; N=%d)', sum(ValidCells4Plot)))
        text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
    end
    subplot(1,2,2)
    scatter(Operant.Info(KDE_OpvsFr_p>0.05),Free.Info(KDE_OpvsFr_p>0.05),40,Colors(KDE_OpvsFr_p>0.05,:), 'filled')
    hold on
    plot([0 7], [0 7], 'r:', 'LineWidth',2)
    hold off
    xlabel('Operant Session')
    ylabel('Free Session')
    ValidCells4Plot = logical((~isnan(Operant.Info .* Free.Info)).*(KDE_OpvsFr_p>0.05));
    if ~Trill
        title(sprintf('Stable Cells with ttest-pvalue>0.05 (bits; N=%d)', sum(ValidCells4Plot)))
        text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
        suplabel(sprintf('Information on %s coherence with sound %s of all calls', MS, FeatureName2))
    elseif Trill==1
        title(sprintf('Stable Cells with ttest-pvalue>0.05 (bits; N=%d)',sum(ValidCells4Plot)))
        text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
        suplabel(sprintf('Information on %s coherence with sound %s of Trills only', MS, FeatureName2))
    end
    
    
    
    figure(8)
    clf
    [h,pstable]=ttest(Operant.Info(KDE_OpvsFr_p>0.05), Free.Info(KDE_OpvsFr_p>0.05))
    [h,pall]=ttest(Operant.Info, Free.Info)
    % [h,pstable]=kstest(Operant.Info(KDE_OpvsFr_p>0.05)-Free.Info(KDE_OpvsFr_p>0.05))
    % [h,pall]=kstest(Operant.Info-Free.Info)
    subplot(1,2,1)
    histogram(Operant.Info - Free.Info,20, 'FaceColor','k')
    xlabel('Information on Motor Coherence difference Operant - Free')
    ylabel('All cells')
    xlim([-5 5])
    hold on
    VL = vline(0, 'r--');
    VL.LineWidth=2;
    title(sprintf('Information change t-test p=%.2e N=%d', pall, sum(~isnan(Operant.Info .* Free.Info))))
    subplot(1,2,2)
    histogram(Operant.Info(KDE_OpvsFr_p>0.05) - Free.Info(KDE_OpvsFr_p>0.05),20, 'FaceColor','k')
    xlim([-5 5])
    hold on
    VL = vline(0, 'r--');
    VL.LineWidth=2;
    xlabel('Information on Motor Coherence difference Operant - Free')
    ylabel('Stable cells')
    title(sprintf('Information change t-test p=%.2e N=%d', pstable, sum((~isnan(Operant.Info .* Free.Info)).*(KDE_OpvsFr_p>0.05))))
    if Trill
        suplabel('Information on coherence of sound Amp with spike rate for Trills only')
    else
        suplabel('Information on coherence of sound Amp with spike rate for all calls only')
    end
    
    
    
    if ~Trill
        MS= 'Auditory';
        Free = load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, 'Free')));
        Operant = load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, 'Operant')));
        Colors = [(Operant.Info_p<0.01) zeros(length(Operant.Info),1) (Free.Info_p<0.01)];
        if strcmp(FeatureName{fn}, 'amp')
            FeatureName2 = 'Amplitude';
        elseif strcmp(FeatureName{fn}, 'sal')
            FeatureName2 = 'Pitch Saliency';
        elseif strcmp(FeatureName{fn}, 'SpectralMean')
            FeatureName2 = 'Spectral Mean';
        end
        figure(6)
        scatter(Operant.Info,Free.Info,40,Colors, 'filled')
        hold on
        plot([0 5], [0 5], 'r:', 'LineWidth',2)
        hold off
        xlabel('Operant Session')
        ylabel('Free Session')
        ValidCells4Plot = ~isnan(Operant.Info .* Free.Info);
        title(sprintf('Information on %s Coherence with sound %s (bits; N=%d)', MS, FeatureName2,sum(~isnan(Operant.Info .* Free.Info))))
        text(0.2, 4.5, sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Operant.Info_p(ValidCells4Plot)<0.01), sum(Free.Info_p(ValidCells4Plot)<0.01),sum(sum(Colors(ValidCells4Plot,:),2)==2)))
    end
% end



%% Compare Information on Coherence between Auditory and Motor in Free session
FeatureName = {'amp' 'SpectralMean' 'sal'};
% for fn = 1:length(FeatureName)
fn=1;

figure(7)
clf
% for fn = 1:length(FeatureName)
 SessionChoice = 'Free';
    Auditory = load(fullfile(Path,sprintf('AuditoryCoherence_%s_%s.mat', FeatureName{fn}, SessionChoice)));
    Motor = load(fullfile(Path,sprintf('MotorCoherence_%s_%s.mat', FeatureName{fn}, SessionChoice)));
    Colors = [(Motor.Info_p<0.01) zeros(length(Info),1) (Auditory.Info_p<0.01)];
    if strcmp(FeatureName{fn}, 'amp')
        FeatureName2 = 'Amplitude';
    elseif strcmp(FeatureName{fn}, 'sal')
        FeatureName2 = 'Pitch Saliency';
    elseif strcmp(FeatureName{fn}, 'SpectralMean')
        FeatureName2 = 'Spectral Mean';
    end
    subplot(2,1,1)
    scatter(Motor.Info,Auditory.Info,40,Colors, 'filled')
    hold on
    plot([0 3], [0 3], 'r:', 'LineWidth',2)
    hold off
    xlabel('Motor')
    ylabel('Auditory')
    ValidCells4Plot = ~isnan(Motor.Info .* Auditory.Info);
    title(sprintf('Information on Coherence with sound %s during %s Session (bits; N=%d)', FeatureName2, SessionChoice, sum(ValidCells4Plot)))
    text(0.2, 2.5,  sprintf('Sig: Red=Motor (N=%d) Blue=Auditory (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
    
    subplot(2,1,2)
    [h,pMA]=ttest(Motor.Info, Auditory.Info)
    histogram(Motor.Info - Auditory.Info,20, 'FaceColor','k')
    xlabel('Information on Coherence difference Motor - Auditory')
    ylabel('All cells')
    xlim([-1.5 1.5])
    hold on
    VL = vline(0, 'r--');
    VL.LineWidth=2;
    title(sprintf('Information change in Free session t-test p=%.2e N=%d', pMA, sum(~isnan(Motor.Info .* Auditory.Info))))
    
% end

%% Effect of sampling on value of coherence
% Cell with Max value of coherence for operant session:
Operant = load(fullfile(Path,sprintf('%sCoherence_%s_%s_Trill.mat', 'Motor', 'amp', 'Operant')));
% CellMaxInd = find(Operant.Info == max(Operant.Info));
CellMaxInd = find((Operant.Info >2.2) .* (Operant.Info <2.4));
CellMaxInd = CellMaxInd(randi(length(CellMaxInd),1));
CellPath = fullfile(Operant.CellsPath(CellMaxInd).folder,CellsPath(CellMaxInd).name);
Cell = load(CellPath, 'What', 'ExpType', 'Who', 'BioSound', 'Duration','SpikesArrivalTimes_Behav');
% find the max number of vocalizations
MaxTrills = length(intersect(intersect(find(contains(Cell.What, 'VocTr')),find(contains(Cell.ExpType, 'O'))), find(contains(Cell.Who, 'self'))));
NumsVoc = 10:5:MaxTrills;
Nrep = 10;
NumVoc_CoherenceInfo = nan(length(NumsVoc),Nrep);
NumVoc_CoherenceInfop = nan(length(NumsVoc),Nrep);
for rr=1:Nrep
    for nv = 1:length(NumsVoc)
        NumVoc = NumsVoc(nv);
        fprintf(1, 'Coherence %d/%d Repetition %d with %d voc\n', nv, length(NumsVoc),rr,NumVoc)
        Coherence_local = coherence4LMC(Cell,1,'Operant', 'amp', 1,NumVoc);
        NumVoc_CoherenceInfo(nv,rr) = Coherence_local.Info;
        NumVoc_CoherenceInfop(nv,rr) = Coherence_local.Info_p;
    end
end

figure()
shadedErrorBar(NumsVoc, mean(NumVoc_CoherenceInfo,2),std(NumVoc_CoherenceInfo,0,2))
hold on
for rr=1:Nrep
    scatter(NumsVoc, NumVoc_CoherenceInfo(:,rr), 40, [NumVoc_CoherenceInfop(:,rr)<0.01 zeros(length(NumVoc_CoherenceInfop(:,rr)),2)], 'filled')
    hold on
end
xlabel('# vocalizations')
ylabel('Information on motor amplitude coherence (bits)')
title(sprintf('Cell: %s', CellPath))
ylim([0 max(max(NumVoc_CoherenceInfo))+0.5])
xlim([0 MaxTrills+5])
yyaxis right
plot(NumsVoc, sum(NumVoc_CoherenceInfop<0.01,2)/Nrep, 'b-')
ylabel('Probability of significance')
ylim([0 1])

hold off

%% Compare Information on Coherence between Rand and True Free session

FeatureName = {'amp' 'SpectralMean' 'sal'};
figure(5)
clf
for fn = 1:length(FeatureName)
    Rand = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s_Rand.mat', FeatureName{fn}, 'Free')));
    Free = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s.mat', FeatureName{fn}, 'Free')));
    if strcmp(FeatureName{fn}, 'amp')
        FeatureName2 = 'Amplitude';
    elseif strcmp(FeatureName{fn}, 'sal')
        FeatureName2 = 'Pitch Saliency';
    elseif strcmp(FeatureName{fn}, 'SpectralMean')
        FeatureName2 = 'Spectral Mean';
    end
    subplot(1,3,fn)
    scatter(Free.Info,Rand.Info,40,[0 0 0], 'filled')
    hold on
    plot([0 5], [0 5], 'r:', 'LineWidth',2)
    hold off
    xlabel('Free Session')
    ylabel('Stim Permutation')
    title(sprintf('Information on Coherence with sound %s (bits)', FeatureName2))
    
    
end

%% Compare Information on Coherence calculated with Trills only between Random and True Operant session

FeatureName = {'amp' 'SpectralMean' 'sal'};
% figure(5)
figure()
clf
% for fn = 1:length(FeatureName)
fn=1;
    Rand = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s_RAND_TRILL.mat', FeatureName{fn}, 'Operant')));
    Operant = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s_TRILL.mat', FeatureName{fn}, 'Operant')));
    if strcmp(FeatureName{fn}, 'amp')
        FeatureName2 = 'Amplitude';
    elseif strcmp(FeatureName{fn}, 'sal')
        FeatureName2 = 'Pitch Saliency';
    elseif strcmp(FeatureName{fn}, 'SpectralMean')
        FeatureName2 = 'Spectral Mean';
    end
%     subplot(1,3,fn)
    scatter(Operant.Info,Rand.Info,40,[0 0 0], 'filled')
    hold on
    plot([0 5], [0 5], 'r:', 'LineWidth',2)
    hold off
    xlabel('Operant Session')
    ylabel('Stim Permutation')
    title(sprintf('Information on Coherence with sound %s using TRILLS only (bits)', FeatureName2))
    
    
% end

%% Compare Information on Coherence calculated with Trills only between Random and True Free session

FeatureName = {'amp' 'SpectralMean' 'sal'};
figure(5)
clf
% for fn = 1:length(FeatureName)
fn=1;
    Rand = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s_RAND_TRILL.mat', FeatureName{fn}, 'Free')));
    Free = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s_TRILL.mat', FeatureName{fn}, 'Free')));
    if strcmp(FeatureName{fn}, 'amp')
        FeatureName2 = 'Amplitude';
    elseif strcmp(FeatureName{fn}, 'sal')
        FeatureName2 = 'Pitch Saliency';
    elseif strcmp(FeatureName{fn}, 'SpectralMean')
        FeatureName2 = 'Spectral Mean';
    end
%     subplot(1,3,fn)
    scatter(Free.Info,Rand.Info,40,[0 0 0], 'filled')
    hold on
    plot([0 5], [0 5], 'r:', 'LineWidth',2)
    hold off
    xlabel('Free Session')
    ylabel('Stim Permutation')
    title(sprintf('Information on Coherence with sound %s using TRILLS only (bits)', FeatureName2))
    
    
% end

%% Compare Information on Coherence calculated with Trills only between Operant and True Free session

FeatureName = {'amp' 'SpectralMean' 'sal'};
figure(5)
clf
% for fn = 1:length(FeatureName)
fn=1;
    Operant = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s_TRILL.mat', FeatureName{fn}, 'Operant')));
    Free = load(fullfile(Path,sprintf('MotorModelsCoherency_%s_%s_TRILL.mat', FeatureName{fn}, 'Free')));
    if strcmp(FeatureName{fn}, 'amp')
        FeatureName2 = 'Amplitude';
    elseif strcmp(FeatureName{fn}, 'sal')
        FeatureName2 = 'Pitch Saliency';
    elseif strcmp(FeatureName{fn}, 'SpectralMean')
        FeatureName2 = 'Spectral Mean';
    end
%     subplot(1,3,fn)
    scatter(Operant.Info,Free.Info,40,[0 0 0], 'filled')
    hold on
    plot([0 7], [0 7], 'r:', 'LineWidth',2)
    hold off
    xlabel('Operant Session')
    ylabel('Free session')
    title(sprintf('Information on Coherence with sound %s using TRILLS only (bits)', FeatureName2))
    xlim([0 7])
    ylim([0 7])
    
% end

%% Explore Cell by cell the profile of coherence and scatter plots of cell tuning for Good Cells
%% Then for all cells
InputFig=0;
load(fullfile(Path,'MotorModelsCoherency_amp.mat'))
GoodInfo=find(~isnan(Info));
NCells = length(GoodInfo);
AmpPredictor = nan(NCells,3); % first element is the F statistic with full model, second is pvalue, third is the partial adjusted R2
SalPredictor = nan(NCells,3);% first element is the F statistic with full model, second is pvalue, third is the partial adjusted R2
SpecMeanPredictor = nan(NCells,3); % first element is the F statistic with full model, second is pvalue, third is the partial adjusted R2
CallTypePredictor = nan(NCells,3); % first element is the F statistic with full model, second is pvalue, third is the partial adjusted R2
FullModelR2 = nan(NCells,3); % first element is the F statistic with Null model, second is pvalue, third is the partial adjusted R2
ModelAmpR2 = nan(NCells,3); % first element is the F statistic with Null model, second is pvalue, third is the adjusted R2, Fourth is the coefficient in the model
ModelSalR2 = nan(NCells,3); % first element is the F statistic with Null model, second is pvalue, third is the adjusted R2, Fourth is the coefficient in the model
ModelSpecMeanR2 = nan(NCells,3); % first element is the F statistic with Null model, second is pvalue, third is the adjusted R2, Fourth is the coefficient in the model
ModelCallTypeR2 = nan(NCells,3); % first element is the F statistic with Null model, second is pvalue, third is the adjusted R2, Fourth is the coefficient in the model
for nc=1:NCells
    cc=GoodInfo(nc);
    fprintf(1,'Cell %d/%d\n',nc,NCells)
    % Figure on coherence values with sound Amplitude
    figure(3)
    clf
    subplot(1,3,1)
    plot(CoherencyT_xTimeDelay{cc}, CoherencyT_filt{cc}, 'LineWidth',2)
    if ~isnan(CoherencyT_DelayAtzero(cc))
        hold on
        P = CoherencyT_filt{cc}(CoherencyT_xTimeDelay{cc} == CoherencyT_DelayAtzero(cc));
        plot(CoherencyT_DelayAtzero(cc),P,'ro', 'MarkerSize',16, 'MarkerFaceColor','r')
        text(max(CoherencyT_xTimeDelay{cc})/5, 4.5/5*P, sprintf('Delay at zero = %.2f ms', CoherencyT_DelayAtzero(cc)))
        text(max(CoherencyT_xTimeDelay{cc})/5, 4/5*P, sprintf('Width at y=mean(coherencyT) = %.2f ms', CoherencyT_WidthAtMaxPeak(cc)))
        hold off
    end
    xlabel('Time Delay')
    ylabel('Coherency')
    
    subplot(1,3,2)
    shadedErrorBar(Freqs{cc}, Coherence{cc},[(Coherence_up{cc}-Coherence{cc})'; (-Coherence_low{cc}+Coherence{cc})'], {'LineWidth',2,'Color','k'})
    hold on
    plot([Freqs{cc}(1) Freqs{cc}(end)], [0 0], 'r--', 'LineWidth',2)
    hold on
    plot(CoherencePeaksF{cc}, CoherencePeaks{cc}, 'go', 'MarkerSize',6,'MarkerFaceColor','g')
    hold on
    plot(MaxCoherence(cc,2), MaxCoherence(cc,1), 'ro', 'MarkerSize',6,'MarkerFaceColor','r')
    MaxY = 0.5;
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
        suplabel(sprintf('Cell %s', CellsPath(cc).name),'t');
    end
    
    % Now plot the scatter tuning curves
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
            fprintf(1,'*** . Problem with Cell %d, duration inconcistency!! ****\n', cc)
            continue
        end
    end
    
    % Get the type of the vocalization
    Trill1_Bark0 = contains(Cell.What(IndVoc), 'Tr');
    
    % Get the optimal Time resolution values
    TRs=[];
    TRs(1) = CoherencyT_WidthAtMaxPeak(cc);
    %     if ~isnan(SecondCoherenceFreqCutOff(cc))
    %         TRs(2) = round(1/SecondCoherenceFreqCutOff(cc)*10^3);
    %     end
    
    
    for tr = 1:length(TRs)
        TR = TRs(tr);
        % Compute neural vectors
        % Neural Data loop
        % neural response is a vector that compile all spike counts starting
        % at 200ms (-Delay(1)) before stim onset and stop at 200ms (Delay(2)) after
        % stim offset
        Delay = [-CoherencyT_DelayAtzero(cc) CoherencyT_DelayAtzero(cc)];
        [YPerStim, YPerStimt] = get_y_4Coherence(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),Delay,TR, TR/2);
        Y = [YPerStim{:}]';
        if isempty(Y)
            fprintf(1,'no spike during vocalization! No model!\n')
            continue
        end
        % Get the vector of call type
        Ind = [0 cumsum(cellfun('length',YPerStim))];
        Trill1_Bark0_local = nan(size(Y));
        for ss = 1:NStims
            Trill1_Bark0_local((Ind(ss)+1):Ind(ss+1)) = Trill1_Bark0(ss) .* ones(Ind(ss+1)-Ind(ss),1);
        end
        
        
        % Calculate acoustic features input to the models
        % acoustic data is a vector of the value of the acoustic feature sampled
        % at 1000Hz starting 200ms (-Delay(1)) before stim onset and stop at 200ms (Delay(2)) after
        % stim offset
        Delay = [0 0];
        DefaultVal = 0;%zero should be the default value for the amplitude, we know here that there is no sound
        [XAmpPerStim,XAmpPerStimt]  = get_x_4coherence(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Delay,TR,DefaultVal,'amp', TR/2);
        XAmp = [XAmpPerStim{:}]';
        
        DefaultVal = 0;%zero should be the default value for the saliency, we know here that there is no sound
        [XSalPerStim, XSalPerStimt] = get_x_4coherence(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Delay,TR,DefaultVal,'sal', TR/2);
        XSal = [XSalPerStim{:}]';
        
        DefaultVal = 'mean';%mean of specmean should be the default value for the saliency, we know here that there is no sound
        [XSpecMeanPerStim, XSpecMeanPerStimt] = get_x_4coherence(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Delay,TR,DefaultVal,'SpectralMean', TR/2);
        XSpecMean = [XSpecMeanPerStim{:}]';
        XSpecMean(isnan(XSpecMean)) = nanmean(XSpecMean);
        
        if InputFig
            plotxyfeaturescoherence(Cell.BioSound(IndVoc,1),YPerStim,YPerStimt,XSpecMeanPerStim,XSpecMeanPerStimt,TR,Delay,Cell.Duration(IndVoc), 50000,'Spectral Mean')
            plotxyfeaturescoherence(Cell.BioSound(IndVoc,2),YPerStim,YPerStimt,XSalPerStim,XSalPerStimt,TR,Delay,Cell.Duration(IndVoc), 10000,'Saliency')
            plotxyfeaturescoherence(Cell.BioSound(IndVoc,2),YPerStim,YPerStimt,XAmpPerStim,XAmpPerStimt,TR,Delay,Cell.Duration(IndVoc), 10000,'Amplitude')
        end
        
        % Run some GLM to establish the effect of acoustic features
        % formula Y~ XAmp + XSal + XSpecMean
        %         TermsMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 0 ];
        TableFull = table(XAmp, XSal, XSpecMean, categorical(Trill1_Bark0_local),Y,'VariableNames',{'Amplitude' 'Saliency' 'SpectralMean' 'CallType' 'Rate'});
        FullModel = fitlm(TableFull, 'Rate ~ Amplitude + Saliency + SpectralMean + CallType + Amplitude:Saliency + Amplitude:SpectralMean + Amplitude:CallType + Saliency:SpectralMean + Saliency:CallType + SpectralMean:CallType' , 'CategoricalVars', 4);
        FullModelR2(nc,3) = FullModel.Rsquared.Adjusted;
        AnovaF = anova(FullModel, 'summary');
        FullModelR2(nc,1) = AnovaF.F(contains(AnovaF.Properties.RowNames, 'Model'));
        FullModelR2(nc,2) = AnovaF.pValue(contains(AnovaF.Properties.RowNames, 'Model'));
        
        ModelwoAmp = fitlm(TableFull, 'Rate ~ Saliency + SpectralMean + CallType + Saliency:SpectralMean + Saliency:CallType + SpectralMean:CallType' , 'CategoricalVars', 4);
        ModelAmp = fitlm(XAmp, Y );
        ModelAmpR2(nc,3) = ModelAmp.Rsquared.Adjusted;
        AnovaAmp = anova(ModelAmp, 'summary');
        ModelAmpR2(nc,1) = AnovaAmp.F(contains(AnovaAmp.Properties.RowNames, 'Model'));
        ModelAmpR2(nc,2) = AnovaAmp.pValue(contains(AnovaAmp.Properties.RowNames, 'Model'));
        if ModelAmp.Rsquared.Adjusted>0
            ModelAmpR2(nc,4) =  ModelAmp.Coefficients.Estimate(2);
        end
        XPredictAmp = min(XAmp):(max(XAmp)-min(XAmp))/10:max(XAmp);
        [YPredictAmp,YPredictAmpci] = predict(ModelAmp,XPredictAmp');
        
        ModelwoSal = fitlm(TableFull, 'Rate ~ Amplitude + SpectralMean + CallType +  Amplitude:SpectralMean + Amplitude:CallType +  SpectralMean:CallType' , 'CategoricalVars', 4);
        ModelSal = fitlm(XSal,Y);
        ModelSalR2(nc,3) = ModelSal.Rsquared.Adjusted;
        AnovaSal = anova(ModelSal, 'summary');
        ModelSalR2(nc,1) = AnovaSal.F(contains(AnovaSal.Properties.RowNames, 'Model'));
        ModelSalR2(nc,2) = AnovaSal.pValue(contains(AnovaSal.Properties.RowNames, 'Model'));
        if ModelSal.Rsquared.Adjusted>0
            ModelSalR2(nc,4) =  ModelSal.Coefficients.Estimate(2);
        end
        XPredictSal = min(XSal):(max(XSal)-min(XSal))/10:max(XSal);
        [YPredictSal,YPredictSalci] = predict(ModelSal,XPredictSal');
        
        ModelwoSpecMean = fitlm(TableFull, 'Rate ~ Amplitude + Saliency + CallType + Amplitude:Saliency +  Amplitude:CallType +  Saliency:CallType' , 'CategoricalVars', 4);
        ModelSpecMean = fitlm(XSpecMean,Y);
        ModelSpecMeanR2(nc,3) = ModelSpecMean.Rsquared.Adjusted;
        AnovaSpecMean = anova(ModelSpecMean, 'summary');
        ModelSpecMeanR2(nc,1) = AnovaSpecMean.F(contains(AnovaSpecMean.Properties.RowNames, 'Model'));
        ModelSpecMeanR2(nc,2) = AnovaSpecMean.pValue(contains(AnovaSpecMean.Properties.RowNames, 'Model'));
        if ModelSpecMean.Rsquared.Adjusted>0
            ModelSpecMeanR2(nc,4) =  ModelSpecMean.Coefficients.Estimate(2);
        end
        XPredictSpecMean = min(XSpecMean):(max(XSpecMean)-min(XSpecMean))/10:max(XSpecMean);
        [YPredictSpecMean,YPredictSpecMeanci] = predict(ModelSpecMean,XPredictSpecMean');
        
        ModelwoCallType = fitlm(TableFull, 'Rate ~ Amplitude + Saliency + SpectralMean +  Amplitude:Saliency + Amplitude:SpectralMean + Saliency:SpectralMean ' , 'CategoricalVars', 4);
        ModelCallType = fitlm(categorical(Trill1_Bark0_local),Y);
        ModelCallTypeR2(nc,3) = ModelCallType.Rsquared.Adjusted;
        AnovaCallType = anova(ModelCallType, 'summary');
        ModelCallTypeR2(nc,1) = AnovaCallType.F(contains(AnovaCallType.Properties.RowNames, 'Model'));
        ModelCallTypeR2(nc,2) = AnovaCallType.pValue(contains(AnovaCallType.Properties.RowNames, 'Model'));
        if ModelCallType.Rsquared.Adjusted>0
            ModelCallTypeR2(nc,4) =  ModelCallType.Coefficients.Estimate(2);
        end
        
        % calculate the significance of the difference in Error according
        % to F distribution and the adjusted R2
        if FullModel.NumObservations ~= ModelwoAmp.NumObservations
            fprintf('Different number of observations')
            keyboard
        end
        kAmp = ModelwoAmp.NumEstimatedCoefficients;
        kFull = FullModel.NumEstimatedCoefficients;
        AmpPredictor(nc,1) = ((ModelwoAmp.SSE - FullModel.SSE)/(kFull - kAmp))/(FullModel.SSE/(FullModel.NumObservations - kFull));
        AmpPredictor(nc,2) = fcdf(AmpPredictor(nc,1), kFull - kAmp, (FullModel.NumObservations - kFull), 'upper');
        AmpPredictor(nc,3) = 1 - (FullModel.SSE / (FullModel.NumObservations - FullModel.NumEstimatedCoefficients))/(ModelwoAmp.SSE / (ModelwoAmp.NumObservations - ModelwoAmp.NumEstimatedCoefficients));
        
        if FullModel.NumObservations ~= ModelwoSal.NumObservations
            fprintf('Different number of observations')
            keyboard
        end
        kSal = ModelwoSal.NumEstimatedCoefficients;
        SalPredictor(nc,1) = ((ModelwoSal.SSE - FullModel.SSE)/(kFull - kSal))/(FullModel.SSE/(FullModel.NumObservations - kFull));
        SalPredictor(nc,2) = fcdf(SalPredictor(nc,1), kFull - kSal, (FullModel.NumObservations - kFull), 'upper');
        SalPredictor(nc,3) = 1 - (FullModel.SSE / (FullModel.NumObservations - FullModel.NumEstimatedCoefficients))/(ModelwoSal.SSE / (ModelwoSal.NumObservations - ModelwoSal.NumEstimatedCoefficients));
        
        if FullModel.NumObservations ~= ModelwoSpecMean.NumObservations
            fprintf('Different number of observations')
            keyboard
        end
        kSpecMean = ModelwoSpecMean.NumEstimatedCoefficients;
        SpecMeanPredictor(nc,1) = ((ModelwoSpecMean.SSE - FullModel.SSE)/(kFull - kSpecMean))/(FullModel.SSE/(FullModel.NumObservations - kFull));
        SpecMeanPredictor(nc,2) = fcdf(SpecMeanPredictor(nc,1), kFull - kSpecMean, (FullModel.NumObservations - kFull), 'upper');
        SpecMeanPredictor(nc,3) = 1 - (FullModel.SSE / (FullModel.NumObservations - FullModel.NumEstimatedCoefficients))/(ModelwoSpecMean.SSE / (ModelwoSpecMean.NumObservations - ModelwoSpecMean.NumEstimatedCoefficients));
        
        if FullModel.NumObservations ~= ModelwoCallType.NumObservations
            fprintf('Different number of observations')
            keyboard
        end
        kCallType = ModelwoCallType.NumEstimatedCoefficients;
        CallTypePredictor(nc,1) = ((ModelwoCallType.SSE - FullModel.SSE)/(kFull - kCallType))/(FullModel.SSE/(FullModel.NumObservations - kFull));
        CallTypePredictor(nc,2) = fcdf(CallTypePredictor(nc,1), kFull - kCallType, (FullModel.NumObservations - kFull), 'upper');
        CallTypePredictor(nc,3) = 1 - (FullModel.SSE / (FullModel.NumObservations - FullModel.NumEstimatedCoefficients))/(ModelwoCallType.SSE / (ModelwoCallType.NumObservations - ModelwoCallType.NumEstimatedCoefficients));
        
        figure(3+tr)
        clf
        subplot(1,3,1)
        shadedErrorBar(XPredictAmp,YPredictAmp/(10^-3),[YPredictAmp-YPredictAmpci(:,1) YPredictAmpci(:,2)- YPredictAmp]'./(10^-3))
        hold on
        scatter(XAmp, Y./(10^-3), 20,[Trill1_Bark0_local zeros(size(Trill1_Bark0_local)) zeros(size(Trill1_Bark0_local))],'filled')
        YLim1 = get(gca, 'YLim');
        XLim1 = get(gca, 'XLim');
        text(XLim1(2)*3/4, YLim1(2)*9.5/10,'Trill','Color',[1 0 0], 'FontWeight','bold')
        text(XLim1(2)*3/4, YLim1(2)*9/10,'Bark','Color',[0 0 0], 'FontWeight','bold')
        xlabel('Sound Amplitude')
        ylabel('Spike Rate (Hz)')
        title(sprintf('R2 = %.2f Partial R2 = %.2f p=%.2f', ModelAmpR2(nc,3), AmpPredictor(nc,[3 2])))
        hold off
        
        subplot(1,3,2)
        shadedErrorBar(XPredictSal,YPredictSal/(10^-3),[YPredictSal-YPredictSalci(:,1) YPredictSalci(:,2)- YPredictSal]'./(10^-3))
        hold on
        scatter(XSal, Y./(10^-3), 20,[Trill1_Bark0_local zeros(size(Trill1_Bark0_local)) zeros(size(Trill1_Bark0_local))],'filled')
        YLim2 = get(gca, 'YLim');
        XLim2 = get(gca, 'XLim');
        text(XLim2(2)*3/4, YLim2(2)*9.5/10,'Trill','Color',[1 0 0], 'FontWeight','bold')
        text(XLim2(2)*3/4, YLim2(2)*9/10,'Bark','Color',[0 0 0], 'FontWeight','bold')
        xlabel('Sound Pitch Saliency')
        ylabel('Spike Rate (Hz)')
        title(sprintf('R2 = %.2f Partial R2 = %.2f p=%.2f', ModelSalR2(nc,3),SalPredictor(nc,[3 2])))
        hold off
        
        subplot(1,3,3)
        shadedErrorBar(XPredictSpecMean,YPredictSpecMean/(10^-3),[YPredictSpecMean-YPredictSpecMeanci(:,1) YPredictSpecMeanci(:,2)- YPredictSpecMean]'./(10^-3))
        hold on
        scatter(XSpecMean, Y./(10^-3), 20,[Trill1_Bark0_local zeros(size(Trill1_Bark0_local)) zeros(size(Trill1_Bark0_local))],'filled')
        YLim3 = get(gca, 'YLim');
        XLim3 = get(gca, 'XLim');
        text(diff(XLim3)*3/4 + XLim3(1), YLim3(2)*9.5/10,'Trill','Color',[1 0 0], 'FontWeight','bold')
        text(diff(XLim3)*3/4 + XLim3(1), YLim3(2)*9/10,'Bark','Color',[0 0 0], 'FontWeight','bold')
        xlabel('Sound Spectral Mean')
        ylabel('Spike Rate (Hz)')
        title(sprintf('R2 = %.2f Partial R2 = %.2f p=%.2f', ModelSpecMeanR2(nc,3), SpecMeanPredictor(nc,[3 2])))
        hold off
        
        suplabel(sprintf('Cell %s TR = %d  ms R2 = %.2f', CellsPath(cc).name, TR,FullModelR2(nc,3)),'t');
        suplabel(sprintf('CallType R2 = %.2f Partial R2 = %.2f p=%.2f', ModelCallTypeR2(nc,3), CallTypePredictor(nc,[3 2])),'x');
        
    end
    pause(1)
    
end

save(fullfile(Path, 'LM_Acoustic.mat'),'CellsPath','GoodInfo','NCells', 'AmpPredictor','SalPredictor','SpecMeanPredictor' , 'CallTypePredictor','FullModelR2','ModelAmpR2','ModelSalR2','ModelSpecMeanR2','ModelCallTypeR2');
%% Plot results of models

% Get the vector of single vs multi units
SU1MU0 = nan(size(CellsPath));
for cc=1:length(CellsPath)
    SU1MU0(cc) = contains(CellsPath(cc).name, 'SSS');
end



% first figure of the R2 values for the full model with all cells
figure(31)
clf
subplot(2,2,1)
scatter(Info(GoodInfo), FullModelR2(:,3), 30, [SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')

subplot(2,2,2)
scatter(Info(GoodInfo), FullModelR2(:,3), 30, [FullModelR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Significant','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')

subplot(2,2,3)
histogram(FullModelR2(:,3), 'BinWidth',0.01,'FaceColor','k')
ylabel('# High Info Cells')
xlabel('Adjusted R2 full linear model')
suplabel('Full acoustic model performance','t')

% Now figure of individual models with all cells
figure(32)
clf
subplot(3,4,1)
scatter(FullModelR2(:,3), ModelAmpR2(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Amplitude')

subplot(3,4,2)
scatter(FullModelR2(:,3), ModelSalR2(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Pitch Saliency')

subplot(3,4,3)
scatter(FullModelR2(:,3), ModelSpecMeanR2(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Spectral Mean')

subplot(3,4,4)
scatter(FullModelR2(:,3), ModelCallTypeR2(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 CallType')

subplot(3,4,5)
scatter(FullModelR2(:,3), AmpPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Amplitude')

subplot(3,4,6)
scatter(FullModelR2(:,3), SalPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,7)
scatter(FullModelR2(:,3), SpecMeanPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,8)
scatter(FullModelR2(:,3), CallTypePredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Full Model R2')

subplot(3,4,9)
scatter(ModelAmpR2(:,4), AmpPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 1])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
hold on
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Amplitude')
ylabel('Partial R2 Amplitude')

subplot(3,4,10)
scatter(ModelSalR2(:,4), SalPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 1])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
hold on

text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Pitch Saliency')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,11)
scatter(ModelSpecMeanR2(:,4), SpecMeanPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 1])
hold on
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Spectral Mean')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,12)
scatter(ModelCallTypeR2(:,4), CallTypePredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 1])
hold on
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Coefficient Call Type')

suplabel('Adjusted R-squared','t')


% Same figure as 32 but now color coded is significance
figure(33)
clf
subplot(3,4,1)
scatter(FullModelR2(:,3), ModelAmpR2(:,3),30,[ModelAmpR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Amplitude')

subplot(3,4,2)
scatter(FullModelR2(:,3), ModelSalR2(:,3),30,[ModelSalR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Pitch Saliency')

subplot(3,4,3)
scatter(FullModelR2(:,3), ModelSpecMeanR2(:,3),30,[ModelSpecMeanR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Spectral Mean')

subplot(3,4,4)
scatter(FullModelR2(:,3), ModelCallTypeR2(:,3),30,[ModelCallTypeR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 CallType')

subplot(3,4,5)
scatter(FullModelR2(:,3), AmpPredictor(:,3),30,[AmpPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Amplitude')

subplot(3,4,6)
scatter(FullModelR2(:,3), SalPredictor(:,3),30,[SalPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,7)
scatter(FullModelR2(:,3), SpecMeanPredictor(:,3),30,[SpecMeanPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,8)
scatter(FullModelR2(:,3), CallTypePredictor(:,3),30,[CallTypePredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Full Model R2')

subplot(3,4,9)
scatter(ModelAmpR2(:,4), AmpPredictor(:,3),30,[AmpPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 1])
hold on
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Amplitude')
ylabel('Partial R2 Amplitude')

subplot(3,4,10)
scatter(ModelSalR2(:,4), SalPredictor(:,3),30,[SalPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 1])
hold on
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Pitch Saliency')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,11)
scatter(ModelSpecMeanR2(:,4), SpecMeanPredictor(:,3),30,[SpecMeanPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 1])
hold on
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Spectral Mean')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,12)
scatter(ModelCallTypeR2(:,4), CallTypePredictor(:,3),30,[CallTypePredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 1])
hold on
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Coefficient Call Type')

suplabel('Adjusted R-squared','t')

% Only keep cells that are significant for at least one model (Full or
% individual feature)

SigCells = find(logical((FullModelR2(:,2)<0.01) + (ModelCallTypeR2(:,2)<0.01) + (ModelSalR2(:,2)<0.01) + (ModelSpecMeanR2(:,2)< 0.01) + (ModelAmpR2(:,2)<0.01)));
% first figure of the R2 values for the full model with all significant cells
figure(34)
clf
subplot(2,2,1)
scatter(Info(GoodInfo(SigCells)), FullModelR2((SigCells),3), 60, [SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')

subplot(2,2,2)
scatter(Info(GoodInfo(SigCells)), FullModelR2((SigCells),3), 60, [FullModelR2((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Significant','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')



subplot(2,2,3)
histogram(FullModelR2((SigCells),3), 'BinWidth',0.01,'FaceColor','k')
ylabel('# High Info Cells')
xlabel('Adjusted R2 full linear model')

subplot(2,2,4)
SigCellsCT = SigCells(CallTypePredictor(SigCells,2)<0.01);
NSigCellsCT = SigCells(CallTypePredictor(SigCells,2)>=0.01);
scatter(Info(GoodInfo(SigCellsCT)), FullModelR2((SigCellsCT),3), 60, [0.8*(AmpPredictor((SigCellsCT),2)<0.01) 0.8*(SalPredictor(SigCellsCT,2)<0.01) 0.8*(SpecMeanPredictor(SigCellsCT,2)<0.01)], 'filled','d')
legend('Call-Type')
legend('AutoUpdate', 'off')
hold on
scatter(Info(GoodInfo(NSigCellsCT)), FullModelR2((NSigCellsCT),3), 60, [0.8*(AmpPredictor((NSigCellsCT),2)<0.01) 0.8*(SalPredictor(NSigCellsCT,2)<0.01) 0.8*(SpecMeanPredictor(NSigCellsCT,2)<0.01)], 'filled','o')

xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*8.5/10,'Amp','Color',[0.8 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*8/10,'Sal','Color',[0 0.8 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*7.5/10,'SpecMean','Color',[0 0 0.8], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*7/10,'Amp + Sal','Color',[0.8 0.8 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*6.5/10,'Amp + SpecMean','Color',[0.8 0 0.8], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*6/10,'Sal + SpecMean','Color',[0 0.8 0.8], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*5.5/10,'Amp + Sal + SpecMean','Color',[0.8 0.8 0.8], 'FontWeight','bold')

suplabel('Full acoustic model performance Significant Cells','t')

% same plot of the population with the significance of the various models
% as color codes on the left and the belonging to a cluster from the Trill
% vs Bark analysis on the right
figure(100)
subplot(1,2,1)
SigCellsCT = SigCells(CallTypePredictor(SigCells,2)<0.01);
NSigCellsCT = SigCells(CallTypePredictor(SigCells,2)>=0.01);
scatter(Info(GoodInfo(SigCellsCT)), FullModelR2((SigCellsCT),3), 60, [0.8*(AmpPredictor((SigCellsCT),2)<0.01) 0.8*(SalPredictor(SigCellsCT,2)<0.01) 0.8*(SpecMeanPredictor(SigCellsCT,2)<0.01)], 'filled','d')
legend('Call-Type')
legend('AutoUpdate', 'off')
hold on
scatter(Info(GoodInfo(NSigCellsCT)), FullModelR2((NSigCellsCT),3), 60, [0.8*(AmpPredictor((NSigCellsCT),2)<0.01) 0.8*(SalPredictor(NSigCellsCT,2)<0.01) 0.8*(SpecMeanPredictor(NSigCellsCT,2)<0.01)], 'filled','o')

xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*8.5/10,'Amp','Color',[0.8 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*8/10,'Sal','Color',[0 0.8 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*7.5/10,'SpecMean','Color',[0 0 0.8], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*7/10,'Amp + Sal','Color',[0.8 0.8 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*6.5/10,'Amp + SpecMean','Color',[0.8 0 0.8], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*6/10,'Sal + SpecMean','Color',[0 0.8 0.8], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*5.5/10,'Amp + Sal + SpecMean','Color',[0.8 0.8 0.8], 'FontWeight','bold')

% retrieve the cluster number
ColorCode = get(groot, 'DefaultAxesColorOrder');
% CellsPath(GoodInfo(SigCells))
ClustNum = nan(size(SigCells));
ClustCol = nan(length(SigCells),3);
ClusRes = load('explore_populationSU_data_5000.mat');
UClust = unique(ClusRes.T);
for cc=1:length(SigCells)
    for uu=1:length(UClust)
        Clust = UClust(uu);
        ClustCellIndices = ClusRes.GoodCellIndices(ClusRes.BaTr_ind(ClusRes.T==Clust));
        CellClust = cell(length(ClustCellIndices),1);
        for ii = 1:length(ClustCellIndices)
            CellClust{ii} = ClusRes.ListSSU(ClustCellIndices(ii)).name;
        end
        if sum(contains(CellClust,CellsPath(GoodInfo(SigCells(cc))).name))
            ClustNum(cc) = Clust;
            ClustCol(cc,:) = ColorCode(Clust,:);
            break
        end
    end
end

subplot(1,2,2)
scatter(Info(GoodInfo(SigCells)), FullModelR2((SigCells),3), 60, ClustCol, 'filled','o')
legend('AutoUpdate', 'off')
hold on
xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
for uu=1:length(UClust)
    text(diff(XLim)*3/4 + XLim(1), YLim(2)*(9-0.5*uu)/10,sprintf('Cluster%d', UClust(uu)),'Color',ColorCode(UClust(uu),:), 'FontWeight','bold')
end
hold off

% Now figure of individual models with all significant cells
figure(36)
clf
subplot(3,4,1)
scatter(FullModelR2((SigCells),3), ModelAmpR2((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Amplitude')

subplot(3,4,2)
scatter(FullModelR2((SigCells),3), ModelSalR2((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Pitch Saliency')

subplot(3,4,3)
scatter(FullModelR2((SigCells),3), ModelSpecMeanR2((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Spectral Mean')

subplot(3,4,4)
scatter(FullModelR2((SigCells),3), ModelCallTypeR2((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 CallType')

subplot(3,4,5)
scatter(FullModelR2((SigCells),3), AmpPredictor((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Amplitude')

subplot(3,4,6)
scatter(FullModelR2((SigCells),3), SalPredictor((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,7)
scatter(FullModelR2((SigCells),3), SpecMeanPredictor((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,8)
scatter(FullModelR2((SigCells),3), CallTypePredictor((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Full Model R2')

subplot(3,4,9)
scatter(ModelAmpR2((SigCells),4), AmpPredictor((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
ylim([-0.2 1])
hold on
plot([0 0], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Amplitude')
ylabel('Partial R2 Amplitude')

subplot(3,4,10)
scatter(ModelSalR2((SigCells),4), SalPredictor((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
ylim([-0.2 1])
hold on
plot([0 0], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Pitch Saliency')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,11)
scatter(ModelSpecMeanR2((SigCells),4), SpecMeanPredictor((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
ylim([-0.2 1])
hold on
plot([0 0], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Spectral Mean')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,12)
scatter(ModelCallTypeR2((SigCells),4), CallTypePredictor((SigCells),3),30,[SU1MU0(GoodInfo(SigCells)) zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
ylim([-0.2 1])
hold on
plot([0 0], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Coefficient Call Type')

suplabel('Adjusted R-squared Significant Cells','t')


% Same figure as 36 but now color coded is significance
figure(37)
clf
subplot(3,4,1)
scatter(FullModelR2((SigCells),3), ModelAmpR2((SigCells),3),30,[ModelAmpR2((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Amplitude')

subplot(3,4,2)
scatter(FullModelR2((SigCells),3), ModelSalR2((SigCells),3),30,[ModelSalR2((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Pitch Saliency')

subplot(3,4,3)
scatter(FullModelR2((SigCells),3), ModelSpecMeanR2((SigCells),3),30,[ModelSpecMeanR2((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Spectral Mean')

subplot(3,4,4)
scatter(FullModelR2((SigCells),3), ModelCallTypeR2((SigCells),3),30,[ModelCallTypeR2((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 CallType')

subplot(3,4,5)
scatter(FullModelR2((SigCells),3), AmpPredictor((SigCells),3),30,[AmpPredictor((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Amplitude')

subplot(3,4,6)
scatter(FullModelR2((SigCells),3), SalPredictor((SigCells),3),30,[SalPredictor((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,7)
scatter(FullModelR2((SigCells),3), SpecMeanPredictor((SigCells),3),30,[SpecMeanPredictor((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,8)
scatter(FullModelR2((SigCells),3), CallTypePredictor((SigCells),3),30,[CallTypePredictor((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Full Model R2')

subplot(3,4,9)
scatter(ModelAmpR2((SigCells),4), AmpPredictor((SigCells),3),30,[AmpPredictor((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
ylim([-0.2 1])
hold on
plot([0 0], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Amplitude')
ylabel('Partial R2 Amplitude')

subplot(3,4,10)
scatter(ModelSalR2((SigCells),4), SalPredictor((SigCells),3),30,[SalPredictor((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
ylim([-0.2 1])
hold on
plot([0 0], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Pitch Saliency')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,11)
scatter(ModelSpecMeanR2((SigCells),4), SpecMeanPredictor((SigCells),3),30,[SpecMeanPredictor((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
ylim([-0.2 1])
hold on
plot([0 0], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Spectral Mean')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,12)
scatter(ModelCallTypeR2((SigCells),4), CallTypePredictor((SigCells),3),30,[CallTypePredictor((SigCells),2)<0.01 zeros(size(SigCells)) zeros(size(SigCells))], 'filled')
ylim([-0.2 1])
hold on
plot([0 0], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Coefficient Call Type')


suplabel('Adjusted R-squared Significant Cells','t')











%% Cells with negative values of Full model R2 are actually not different
% than zero
% Negative values of R2 adj conrespond to 0 values. Adjusted turns negative
% as soon as the ratio of SSE/SSTot > (n-k-1)/(n-1) for instance 0.92 for
% n=100 and k=7
FullModelR2_0corr = FullModelR2(:,3);
FullModelR2_0corr(FullModelR2(:,3)<0) = 0;

% Then cells with higher values of adjusted R2 for single acoustic
% parameter compared to Full Model R2 should be set at the value of Full
% model R2, they are higher just by the mechanism of the adjusted
% correction (more parameters on the full  model)
ModelCallTypeR2_corr = ModelCallTypeR2(:,3);
ModelCallTypeR2_corr(ModelCallTypeR2(:,3)>FullModelR2(:,3)) = FullModelR2_0corr(ModelCallTypeR2(:,3)>FullModelR2(:,3));
ModelCallTypeR2_corr(ModelCallTypeR2_corr<0) = 0;

ModelSalR2_corr = ModelSalR2(:,3);
ModelSalR2_corr(ModelSalR2(:,3)>FullModelR2(:,3)) = FullModelR2_0corr(ModelSalR2(:,3)>FullModelR2(:,3));
ModelSalR2_corr(ModelSalR2_corr<0) = 0;

ModelAmpR2_corr = ModelAmpR2(:,3);
ModelAmpR2_corr(ModelAmpR2(:,3)>FullModelR2(:,3)) = FullModelR2_0corr(ModelAmpR2(:,3)>FullModelR2(:,3));
ModelAmpR2_corr(ModelAmpR2_corr<0) = 0;

ModelSpecMeanR2_corr = ModelSpecMeanR2(:,3);
ModelSpecMeanR2_corr(ModelSpecMeanR2(:,3)>FullModelR2(:,3)) = FullModelR2_0corr(ModelSpecMeanR2(:,3)>FullModelR2(:,3));
ModelSpecMeanR2_corr(ModelSpecMeanR2_corr<0) = 0;


% Then cells with higher values of adjusted R2 for restricted Full models by one single acoustic
% parameter compared to Full Model R2 should be set at the value of Full
% model R2, they are higher just by the mechanism of the adjusted
% correction (more parameters on the full  model)
ModelCallTypePartialR2_corr = CallTypePredictor(:,3);
ModelCallTypePartialR2_corr(CallTypePredictor(:,3)>FullModelR2(:,3)) = FullModelR2_0corr(CallTypePredictor(:,3)>FullModelR2(:,3));
ModelCallTypePartialR2_corr(ModelCallTypePartialR2_corr<0) = 0;

ModelSalPartialR2_corr = SalPredictor(:,3);
ModelSalPartialR2_corr(SalPredictor(:,3)>FullModelR2(:,3)) = FullModelR2_0corr(SalPredictor(:,3)>FullModelR2(:,3));
ModelSalPartialR2_corr(ModelSalPartialR2_corr<0) = 0;

ModelAmpPartialR2_corr = AmpPredictor(:,3);
ModelAmpPartialR2_corr(AmpPredictor(:,3)>FullModelR2(:,3)) = FullModelR2_0corr(AmpPredictor(:,3)>FullModelR2(:,3));
ModelAmpPartialR2_corr(ModelAmpPartialR2_corr<0) = 0;

ModelSpecMeanPartialR2_corr = SpecMeanPredictor(:,3);
ModelSpecMeanPartialR2_corr(SpecMeanPredictor(:,3)>FullModelR2(:,3)) = FullModelR2_0corr(SpecMeanPredictor(:,3)>FullModelR2(:,3));
ModelSpecMeanPartialR2_corr(ModelSpecMeanPartialR2_corr<0) = 0;

% Plot again figure 31
figure(34)
clf
subplot(2,2,1)
scatter(Info(GoodInfo), FullModelR2_0corr, 30, [SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')

subplot(2,2,2)
scatter(Info(GoodInfo), FullModelR2_0corr, 30, [FullModelR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Significant','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-significant','Color',[0 0 0], 'FontWeight','bold')

subplot(2,2,3)
histogram(FullModelR2_0corr, 'BinWidth',0.01,'FaceColor','k')
ylabel('# High Info Cells')
xlabel('Adjusted R2 full linear model')
suplabel('Full acoustic model performance','t')

% Plot again figure 32
figure(35)
clf
subplot(3,4,1)
scatter(FullModelR2_0corr, ModelAmpR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Amplitude')

subplot(3,4,2)
scatter(FullModelR2_0corr, ModelSalR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Pitch Saliency')

subplot(3,4,3)
scatter(FullModelR2_0corr, ModelSpecMeanR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Spectral Mean')

subplot(3,4,4)
scatter(FullModelR2_0corr, ModelCallTypeR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 CallType')

subplot(3,4,5)
scatter(FullModelR2_0corr, ModelAmpPartialR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Amplitude')

subplot(3,4,6)
scatter(FullModelR2_0corr, ModelSalPartialR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('R2 Pitch Saliency')

subplot(3,4,7)
scatter(FullModelR2_0corr, ModelSpecMeanPartialR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,8)
scatter(FullModelR2_0corr, ModelCallTypePartialR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Full Model R2')


subplot(3,4,9)
scatter(ModelAmpR2_corr, ModelAmpPartialR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Amplitude')
ylabel('Partial R2 Amplitude')

subplot(3,4,10)
scatter(ModelSalR2_corr, ModelSalPartialR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Saliency')
ylabel('R2 Pitch Saliency')

subplot(3,4,11)
scatter(ModelSpecMeanR2_corr, ModelSpecMeanPartialR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Spectral Mean')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,12)
scatter(ModelCallTypeR2_corr, ModelCallTypePartialR2_corr,30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('R2 Call Type')

suplabel('Corrected Adjusted R-squared','t')


% plot again figure 33
figure(36)
clf
subplot(3,4,1)
scatter(FullModelR2_0corr, ModelAmpR2_corr,30,[ModelAmpR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Amplitude')

subplot(3,4,2)
scatter(FullModelR2_0corr, ModelSalR2_corr,30,[ModelSalR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Pitch Saliency')

subplot(3,4,3)
scatter(FullModelR2_0corr, ModelSpecMeanR2_corr,30,[ModelSpecMeanR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Spectral Mean')

subplot(3,4,4)
scatter(FullModelR2_0corr, ModelCallTypeR2_corr,30,[ModelCallTypeR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 CallType')

subplot(3,4,5)
scatter(FullModelR2_0corr, ModelAmpPartialR2_corr,30,[AmpPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non_Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Amplitude')

subplot(3,4,6)
scatter(FullModelR2_0corr, ModelSalPartialR2_corr,30,[SalPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('R2 Pitch Saliency')

subplot(3,4,7)
scatter(FullModelR2_0corr, ModelSpecMeanPartialR2_corr,30,[SpecMeanPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,8)
scatter(FullModelR2_0corr, ModelCallTypePartialR2_corr,30,[CallTypePredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Full Model R2')


subplot(3,4,9)
scatter(ModelAmpR2_corr, ModelAmpPartialR2_corr,30,[AmpPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Amplitude')
ylabel('Partial R2 Amplitude')

subplot(3,4,10)
scatter(ModelSalR2_corr, ModelSalPartialR2_corr,30,[SalPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Saliency')
ylabel('R2 Pitch Saliency')

subplot(3,4,11)
scatter(ModelSpecMeanR2_corr, ModelSpecMeanPartialR2_corr,30,[SpecMeanPredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Spectral Mean')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,12)
scatter(ModelCallTypeR2_corr, ModelCallTypePartialR2_corr,30,[CallTypePredictor(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 1])
ylim([-0.2 1])
hold on
plot([-0.2 1], [-0.2 1],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'Sig','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'Non-Sig','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('R2 Call Type')

suplabel('Corrected Adjusted R-squared','t')

%% MRFS: Motor Models parameters Ridge regression
load(fullfile(Path,'MotorModelsCoherency_amp.mat'))
% Assumption of stationarity over time

% Running through cells to find the optimal time resolution of the neural response for acoustic feature predicion from the neural response
% here we use a ridge regularization with an MSE optimization on the
% prediction of spike rate
Win =300; %value in ms for the duration of the snippet of sound which
% acoustic features are used to predict the neural response at a given time
% t. This will be the size of the xaxis of the MRF.

Delay = 100;% ms
%The neural activity at time t is predicted by a time varying
% acoustic feature that starts at t-Delay ms
NCells = length(GoodInfo);
% Define the time resolution at which neural density estimates should be calculated
TRs = cell(NCells,1); % time in ms with first column being the optimal time according to coherency peak width and second column being the optimal time according to last peak in coherence
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
    TRs{cc}(1) = CoherencyT_WidthAtMaxPeak(GoodInfo(cc));
    if ~isnan(SecondCoherenceFreqCutOff(GoodInfo(cc)))
        TRs{cc}(2) = round(1/SecondCoherenceFreqCutOff(GoodInfo(cc))*10^3);
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
    MSE_TR_Amp{cc} = nan(1,length(TRs{cc}));
    MSE_TR_SpecMean{cc} = nan(1,length(TRs{cc}));
    MSE_TR_Sal{cc} =  nan(1,length(TRs{cc}));
    MeanYTrain{cc} = nan(1,length(TRs{cc}));
    Ypredict_Amp{cc} = cell(1,length(TRs{cc}));
    Ypredict_SpecMean{cc} = cell(1,length(TRs{cc}));
    Ypredict_Sal{cc} =  cell(1,length(TRs{cc}));
    Yval{cc} =  cell(1,length(TRs{cc}));
    TicToc{cc} = nan(1,length(TRs{cc}));
    Tr1=1;
    
    %%
    for tr = Tr1:length(TRs{cc})
        TR = TRs{cc}(tr);
        fprintf(1,'Cell %d/%d Ridge models with Time resolution %d ms (%d/%d)\n', cc,NCells,TR, tr, length(TRs{cc}));
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
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d) => done in %.2f minutes \n', cc,NCells,TR, tr, length(TRs{cc}), TicToc{cc}(tr)/60);
    end
    if OutFig
        figure(3)
        clf
        plot(TRs{cc}, MSE_TR_Amp{cc}, 'Linewidth',2)
        hold on
        plot(TRs{cc}, MSE_TR_SpecMean{cc}, 'Linewidth',2)
        hold on
        plot(TRs{cc}, MSE_TR_Sal{cc}, 'Linewidth',2)
        xlabel('Time resolution in ms')
        ylabel('Mean Squared Error')
        hold off
        
        figure(4)
        clf
        subplot(2,3,1)
        MAP = colormap();
        c=linspace(1,256,length(TRs{cc}));
        for tr = 1:length(TRs{cc})
            scatter(Yval{cc}, Ypredict_Amp{cc}{tr},30,MAP(round(c(tr)),:),'filled')
            hold on
        end
        xlabel('Observed rate')
        ylabel('Predicted rate')
        title('Amplitude Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs{cc})
        hold off
        YL = ylim;
        XL = xlim;
        ylim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        xlim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        
        subplot(2,3,2)
        for tr = 1:length(TRs{cc})
            scatter(Yval{cc}, Ypredict_SpecMean{cc}{tr},30,MAP(round(c(tr)),:),'filled')
            hold on
        end
        xlabel('Observed rate')
        ylabel('Predicted rate')
        title('SpecMean Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs{cc})
        hold off
        YL = ylim;
        XL = xlim;
        ylim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        xlim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        
        subplot(2,3,3)
        for tr = 1:length(TRs{cc})
            scatter(Yval{cc}, Ypredict_Sal{cc}{tr},30,MAP(round(c(tr)),:),'filled')
            hold on
        end
        xlabel('Observed rate')
        ylabel('Predicted rate')
        title('Saliency Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs{cc})
        hold off
        YL = ylim;
        XL = xlim;
        ylim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        xlim([min(YL(1),XL(1)) max(YL(2),XL(2))])
        
        
        MSE_TR_Amp_local = nan(1,length(TRs{cc}));
        MSE_TR_SpecMean_local = nan(1,length(TRs{cc}));
        MSE_TR_Sal_local = nan(1,length(TRs{cc}));
        subplot(2,3,4)
        MAP = colormap();
        c=linspace(1,256,length(TRs{cc}));
        for tr = 1:length(TRs{cc})
            scatter(Yval{cc}, (Ypredict_Amp{cc}{tr}-Yval{cc}).^2,30,MAP(round(c(tr)),:),'filled')
            hold on
            MSE_TR_Amp_local(tr) = mean((Ypredict_Amp{cc}{tr}-Yval{cc}).^2);
            line(xlim,MSE_TR_Amp_local(tr)*ones(2,1),'Color',MAP(round(c(tr)),:), 'LineWidth',2)
            hold on
        end
        xlabel('Observed rate')
        ylabel('Error2 on rate')
        title('Amplitude Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs{cc})
        hold off
        
        
        subplot(2,3,5)
        for tr = 1:length(TRs{cc})
            scatter(Yval{cc}, (Ypredict_SpecMean{cc}{tr}-Yval{cc}).^2,30,MAP(round(c(tr)),:),'filled')
            hold on
            MSE_TR_SpecMean_local(tr) = mean((Ypredict_SpecMean{cc}{tr}-Yval{cc}).^2);
            line(xlim,MSE_TR_SpecMean_local(tr)*ones(2,1),'Color',MAP(round(c(tr)),:), 'LineWidth',2)
            hold on
        end
        xlabel('Observed rate')
        ylabel('Error2 on log rate')
        title('SpecMean Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs{cc})
        hold off
        
        subplot(2,3,6)
        for tr = 1:length(TRs{cc})
            scatter(Yval{cc}, (Ypredict_Sal{cc}{tr}-Yval{cc}).^2,30,MAP(round(c(tr)),:),'filled')
            hold on
            MSE_TR_Sal_local(tr) = mean((Ypredict_Sal{cc}{tr}-Yval{cc}).^2);
            line(xlim,mean((Ypredict_Sal{cc}{tr}-Yval{cc}).^2)*ones(2,1),'Color',MAP(round(c(tr)),:), 'LineWidth',2)
            hold on
        end
        xlabel('Observed rate')
        ylabel('Error2 on log rate')
        title('Saliency Model')
        colorbar('southoutside','Ticks', c/256,'TickLabels',TRs{cc})
        hold off
        
        figure(3)
        hold on
        plot(TRs{cc}, MSE_TR_Amp_local, 'Linewidth',2, 'LineStyle','--')
        hold on
        plot(TRs{cc}, MSE_TR_SpecMean_local, 'Linewidth',2,'LineStyle','--')
        hold on
        plot(TRs{cc}, MSE_TR_Sal_local, 'Linewidth',2,'LineStyle','--')
        legend({'Amp' 'SpectralMean' 'Saliency' 'AmpMe' 'SpectralMeanMe' 'SaliencyMe'})
        hold off
        
        figure(5)
        clf
        plot(TRs{cc},MeanYTrain{cc}-mean(Yval{cc}))
        xlabel('Time resolution (ms)')
        ylabel('Mean rate difference Training-testing')
    end
    %     keyboard
end
save(fullfile(Path,'MotorModelsRidge.mat'), 'NCells', 'MSE_TR_Amp', 'MSE_TR_SpecMean','MSE_TR_Sal','Ypredict_Amp','Ypredict_SpecMean','Ypredict_Sal','Yval','MeanYTrain','TicToc','CellsPath','TRs', 'GoodInfo');


%% Plot the results of the time resolution optimization using ridge regression
load(fullfile(Path,'MotorModelsRidge.mat'));
% Plot the zscored average MSE accross cells % This does not make sense
% anymore we're already taking best values of Time resolution
% MSE_TR_Amp_zs = nan(NCells,size(TRs,2));
% MSE_TR_SpecMean_zs = nan(NCells,size(TRs,2));
% MSE_TR_Sal_zs = nan(NCells,size(TRs,2));
% also the MSE as a proportion of the average rate of the training correspoding portion
% of the validating set (same temporal resolution)
MSE0 = nan(NCells,2);
R2_TR_Amp = nan(NCells,2);
R2_TR_SpecMean = nan(NCells,2);
R2_TR_Sal = nan(NCells,2);
for cc=1:NCells
    if isempty(MSE_TR_Amp{cc})
        continue
    end
    %     MSE_TR_Amp_zs(cc,:) = zscore(MSE_TR_Amp{cc});
    %     MSE_TR_SpecMean_zs(cc,:) = zscore(MSE_TR_SpecMean{cc});
    %     MSE_TR_Sal_zs(cc,:) = zscore(MSE_TR_Sal{cc});
    for tr=1:length(TRs{cc})
        MSE0(cc,tr) = mean((Yval{cc}-MeanYTrain{cc}(tr)).^2);
    end
    R2_TR_Amp(cc,:) = (MSE0(cc,:)-MSE_TR_Amp{cc})./MSE0(cc,:);
    R2_TR_SpecMean(cc,:) = (MSE0(cc,:)-MSE_TR_SpecMean{cc})./MSE0(cc,:);
    R2_TR_Sal(cc,:) = (MSE0(cc,:)-MSE_TR_Sal{cc})./MSE0(cc,:);
end


% figure(7)
% ColorCode = get(groot, 'DefaultAxesColorOrder');
% subplot(1,2,1)
% legend('AutoUpdate', 'on')
% plot(TRs(cc,:), nanmean(MSE_TR_Amp_zs), 'LineWidth',2, 'Color',ColorCode(1,:), 'DisplayName','Amplitude')
% hold on
% plot(TRs(cc,:), nanmean(MSE_TR_SpecMean_zs), 'LineWidth',2, 'Color',ColorCode(2,:),'DisplayName','SpectralMean')
% hold on
% plot(TRs(cc,:), nanmean(MSE_TR_Sal_zs), 'LineWidth',2, 'Color',ColorCode(3,:),'DisplayName','Saliency')
% legend('AutoUpdate', 'off')
% hold on
% shadedErrorBar(TRs(cc,:), nanmean(MSE_TR_Amp_zs),nanstd(MSE_TR_Amp_zs)./(sum(~isnan(MSE_TR_Amp_zs))).^0.5, {'Color',ColorCode(1,:)})
% hold on
% shadedErrorBar(TRs(cc,:), nanmean(MSE_TR_SpecMean_zs),nanstd(MSE_TR_SpecMean_zs)./(sum(~isnan(MSE_TR_SpecMean_zs))).^0.5, {'Color', ColorCode(2,:)})
% hold on
% shadedErrorBar(TRs(cc,:), nanmean(MSE_TR_Sal_zs),nanstd(MSE_TR_Sal_zs)./(sum(~isnan(MSE_TR_Sal_zs))).^0.5, {'Color', ColorCode(3,:)})
% hold on
% plot(TRs(cc,:), nanmean(MSE_TR_Amp_zs), 'LineWidth',2, 'Color',ColorCode(1,:))
% hold on
% plot(TRs(cc,:), nanmean(MSE_TR_SpecMean_zs), 'LineWidth',2, 'Color',ColorCode(2,:))
% hold on
% plot(TRs(cc,:), nanmean(MSE_TR_Sal_zs), 'LineWidth',2, 'Color',ColorCode(3,:))
% hold off
% xlabel('Time resolution in ms')
% ylabel('zscored Mean Squared Error')
%
% subplot(1,2,2)
% legend('AutoUpdate', 'on')
% plot(TRs(cc,:), nanmean(R2_TR_Amp), 'LineWidth',2, 'Color',ColorCode(1,:), 'DisplayName','Amplitude')
% hold on
% plot(TRs(cc,:), nanmean(R2_TR_SpecMean), 'LineWidth',2, 'Color',ColorCode(2,:),'DisplayName','SpectralMean')
% hold on
% plot(TRs(cc,:), nanmean(R2_TR_Sal), 'LineWidth',2, 'Color',ColorCode(3,:),'DisplayName','Saliency')
% hold on
% legend('AutoUpdate', 'off')
% shadedErrorBar(TRs(cc,:), nanmean(R2_TR_Amp),nanstd(R2_TR_Amp)./(sum(~isnan(R2_TR_Amp))).^0.5, {'Color',ColorCode(1,:)})
% hold on
% shadedErrorBar(TRs(cc,:), nanmean(R2_TR_SpecMean),nanstd(R2_TR_SpecMean)./(sum(~isnan(R2_TR_SpecMean))).^0.5, {'Color', ColorCode(2,:)})
% hold on
% shadedErrorBar(TRs(cc,:), nanmean(R2_TR_Sal),nanstd(R2_TR_Sal)./(sum(~isnan(R2_TR_Sal))).^0.5, {'Color', ColorCode(3,:)})
% hold on
%
% plot(TRs(cc,:), nanmean(R2_TR_Amp), 'LineWidth',2, 'Color',ColorCode(1,:))
% hold on
% plot(TRs(cc,:), nanmean(R2_TR_SpecMean), 'LineWidth',2, 'Color',ColorCode(2,:))
% hold on
% plot(TRs(cc,:), nanmean(R2_TR_Sal), 'LineWidth',2, 'Color',ColorCode(3,:))
% hold off
% xlabel('Time resolution in ms')
% ylabel('R2')
% hold off

% Plot MSE
figure()

ColorCode = get(groot, 'DefaultAxesColorOrder');
subplot(1,3,1)
for cc=1:NCells
    if rem(cc,size(ColorCode,1))
        plot(TRs{cc}, MSE_TR_Amp{cc}, 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(rem(cc,size(ColorCode,1)),:),'MarkerEdgeColor',ColorCode(rem(cc,size(ColorCode,1)),:), 'Color',ColorCode(rem(cc,size(ColorCode,1)),:));
    else
        plot(TRs{cc}, MSE_TR_Amp{cc}, 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(end,:),'MarkerEdgeColor',ColorCode(end,:), 'Color',ColorCode(end,:));
    end
    hold on
end
h=hline(0,'k:');
h.LineWidth = 2;
xlabel('Time Resolution (ms)')
ylabel('MSE Amplitude')


subplot(1,3,2)
for cc=1:NCells
    if rem(cc,size(ColorCode,1))
        plot(TRs{cc}, MSE_TR_SpecMean{cc}, 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(rem(cc,size(ColorCode,1)),:),'MarkerEdgeColor',ColorCode(rem(cc,size(ColorCode,1)),:), 'Color',ColorCode(rem(cc,size(ColorCode,1)),:));
    else
        plot(TRs{cc}, MSE_TR_SpecMean{cc}, 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(end,:),'MarkerEdgeColor',ColorCode(end,:),'Color',ColorCode(end,:));
    end
    hold on
end
h=hline(0,'k:');
h.LineWidth = 2;
xlabel('Time Resolution (ms)')
ylabel('MSE SpectralMean')

subplot(1,3,3)
for cc=1:NCells
    if rem(cc,size(ColorCode,1))
        plot(TRs{cc}, MSE_TR_Sal{cc}, 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(rem(cc,size(ColorCode,1)),:),'MarkerEdgeColor',ColorCode(rem(cc,size(ColorCode,1)),:), 'Color',ColorCode(rem(cc,size(ColorCode,1)),:));
    else
        plot(TRs{cc}, MSE_TR_Sal{cc}, 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(end,:),'MarkerEdgeColor',ColorCode(end,:),'Color',ColorCode(end,:));
    end
    hold on
end
h=hline(0,'k:');
h.LineWidth = 2;
xlabel('Time Resolution (ms)')
ylabel('MSE Saliency')

suplabel('Mean Squared Error Ridge regression models','t')

figure()
ColorCode = get(groot, 'DefaultAxesColorOrder');
subplot(1,3,1)
for cc=1:NCells
    if rem(cc,size(ColorCode,1))
        plot(TRs{cc}, R2_TR_Amp(cc,:), 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(rem(cc,size(ColorCode,1)),:),'MarkerEdgeColor',ColorCode(rem(cc,size(ColorCode,1)),:), 'Color',ColorCode(rem(cc,size(ColorCode,1)),:));
    else
        plot(TRs{cc}, R2_TR_Amp(cc,:), 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(end,:),'MarkerEdgeColor',ColorCode(end,:), 'Color',ColorCode(end,:));
    end
    hold on
end
h=hline(0,'k:');
h.LineWidth = 2;
xlabel('Time Resolution (ms)')
ylabel('R2 Amplitude')


subplot(1,3,2)
for cc=1:NCells
    if rem(cc,size(ColorCode,1))
        plot(TRs{cc}, R2_TR_SpecMean(cc,:), 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(rem(cc,size(ColorCode,1)),:),'MarkerEdgeColor',ColorCode(rem(cc,size(ColorCode,1)),:), 'Color',ColorCode(rem(cc,size(ColorCode,1)),:));
    else
        plot(TRs{cc}, R2_TR_SpecMean(cc,:), 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(end,:),'MarkerEdgeColor',ColorCode(end,:),'Color',ColorCode(end,:));
    end
    hold on
end
h=hline(0,'k:');
h.LineWidth = 2;
xlabel('Time Resolution (ms)')
ylabel('R2 SpectralMean')

subplot(1,3,3)
for cc=1:NCells
    if rem(cc,size(ColorCode,1))
        plot(TRs{cc}, R2_TR_Sal(cc,:), 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(rem(cc,size(ColorCode,1)),:),'MarkerEdgeColor',ColorCode(rem(cc,size(ColorCode,1)),:), 'Color',ColorCode(rem(cc,size(ColorCode,1)),:));
    else
        plot(TRs{cc}, R2_TR_Sal(cc,:), 'o-','MarkerSize',5,'MarkerFaceColor',ColorCode(end,:),'MarkerEdgeColor',ColorCode(end,:),'Color',ColorCode(end,:));
    end
    hold on
end
h=hline(0,'k:');
h.LineWidth = 2;
xlabel('Time Resolution (ms)')
ylabel('R2 Saliency')

suplabel('R2 Ridge regression models','t')




%% MRFS: Run ridge GLM Poisson on acoustic features for cells with high values of Info on coherence (Channel capacity with Amplitude)
% These did not work well...
%load(fullfile(Path,'MotorModelsGLM'),'MotorModels','CellsPath', 'GoodInfo');
load(fullfile(Path,'MotorModelsCoherency_amp.mat'))
% parameters of the Poisson GLM Models
ParamModel.LAMBDARATIO=1e-4;
ParamModel.NUMLAMBDA=10;%25?
ParamModel.LINK='log';
ParamModel.DISTR='poisson';
% Determine a list of alpha (parameter that range the regularization
% betweeen ridge (L2, alpha=0) and lasso (L1, alpha =1))
ParamModel.Alpha=0.001; % STRFs are easier to interpret using ridge than using Lasso and deviances are similar.

NCells = length(GoodInfo);

MotorModels = cell(NCells,1);
parfor cc=1:NCells
    if ~isempty(MotorModels{cc}) % This one was already calculated
        fprintf(1, 'Cell %d/%d Already calculated\n',cc,NCells)
        continue
    end
    Cell = load(fullfile(CellsPath(GoodInfo(cc)).folder,CellsPath(GoodInfo(cc)).name));
    
    
    % Number of vocalizations in the dataset
    if ~isfield(Cell, 'What')
        fprintf(1,'*** . Problem with Cell %d/%d %s, no what field!! ****\n', cc, NCells, CellsPath(GoodInfo(cc)).name)
    end
    IndVoc = find(contains(Cell.What, 'Voc') .* (Cell.Duration>20)); % No saliency calculated when duration <20ms
    NStims = length(IndVoc);
    
    %% Variable organization and running models
    % organize acoustic data as a matrix where each
    % column corresponds to amp env at t-100:Win:t+100, t being
    % time of neural window;
    % neural response is a vector that compile all spike counts starting
    % at -100ms before stim onset and stop at 100ms after
    % stim offset
    BestDevSpecMean = nan(length(TRs{cc}),1); % contains the deviance of the spectral mean model
    BestDevSal = nan(length(TRs{cc}),1);  % contains the deviance of the saliency model
    BestDevAmp = nan(length(TRs{cc}),1);  % contains the deviance of the amplitude model
    % BestDevSpecMed = nan(length(TRs),1); % contains the deviance of the spectral median model
    BestDevNull = nan(length(TRs{cc}),1);  % contains the deviance of the null model
    
    BSpecMean = nan(length(TRs{cc}),2*Win+1); % contains the Betas of the spectral mean model with the best deviance
    BSal = nan(length(TRs{cc}),2*Win+1);  % contains the Betas of the saliency model with the best deviance
    BAmp = nan(length(TRs{cc}),Win+1);  % contains the Betas of the amplitude model with the best deviance
    % BestDevSpecMed = nan(length(TRs),1); % contains the Betas of the spectral
    % median model with the best deviance
    % BNull = nan(length(TRs),Win+1);  % contains the Beta of the null model
    BNull = nan(length(TRs{cc}),1);  % contains the Beta of the null model
    TicToc = nan(length(TRs{cc}),1);  % duration of each loop of models
    MeanY = nan(length(TRs{cc}),1);
    
    for tr = 1:length(TRs{cc})
        TR_local = TRs{cc}(tr);
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d)\n', cc,NCells,TR_local, tr, length(TRs{cc}));
        %% Gather the data
        
        TimerLoopGLM =tic;
        % Neural Data loop
        YPerStim = get_y(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),Win,Delay,TR_local);
        Y_local = [YPerStim{:}]';
        
        MeanY(tr)=nanmean(Y_local);
        % check that the neural response is high enough
        if nanmean(Y_local)<10^-2
            fprintf(1, 'Cell %d/%d Models with Time resolution %d ms (%d/%d) -> Rate is too low, not calculating models to avoid convergence issues\n',cc,NCells,TR_local, tr, length(TRs{cc}))
            continue
        end
        
        
        % Acoustic feature loop, first column of biosound is microphone second is
        % piezo
        XSpecMeanPerStim = get_x(Cell.BioSound(IndVoc,1), Cell.Duration(IndVoc), Win, Delay,'SpectralMean');
        XSpecMean = [XSpecMeanPerStim{:}]';
        XSpecMean(isnan(XSpecMean))=nanmean(reshape(XSpecMean,numel(XSpecMean),1)); % Change the padding with NaN to the mean for the SpecMean for silence.
        XSaliencyPerStim = get_x(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Win, Delay,'sal');
        XSaliency = [XSaliencyPerStim{:}]';
        XSaliency(isnan(XSaliency))=0; % Change the padding with NaN to zeros for the Saliency, we know Silence is not harmonic
        XAmpPerStim = get_x(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Win, Delay,'amp');
        XAmp = [XAmpPerStim{:}]';
        XAmp(isnan(XAmp))=0; % Change the padding with NaN to zeros for the amplitude, we know here that there is no sound
        %         XSpecMedPerStim = get_x(BioSound(IndVoc,1), Duration(IndVoc), Win, TR, Delay,'Q2t');
        %         XSpecMed = [XSpecMedPerStim{:}]';
        
        %% Plot of features and neuronal response if requested (BioSound,YPerStim,XPerStim, TR,Delay,F_high,FeatureName)
        if DatFig
            plotxyfeatures(BioSound(IndVoc,1),YPerStim,XSpecMeanPerStim,TR_local,Win,Delay,Duration(IndVoc), 50000,'Spectral Mean')
            plotxyfeatures(BioSound(IndVoc,2),YPerStim,XSaliencyPerStim,TR_local,Win,Delay,Duration(IndVoc), 10000,'Saliency')
            plotxyfeatures(BioSound(IndVoc,2),YPerStim,XAmpPerStim,TR_local,Win,Delay,Duration(IndVoc), 10000,'Amp')
            %         plotxyfeatures(BioSound(IndVoc,1),YPerStim,XSpecMedPerStim,TR,Delay,Duration(IndVoc), 50000,'Spectral Median')
        end
        
        %% Calculate the model
        % get rid of Nans
        Nan_Ind=find(isnan(Y_local));
        Y_local(Nan_Ind) = [];
        XSpecMean(Nan_Ind,:) =[];
        XSaliency(Nan_Ind,:) =[];
        %         XSpecMed(Nan_Ind,:) =[];
        
        %% Run ridge GLM Poisson on acoustic features
        
        % spectral mean predicting Y
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d): SpectralMean Model\n', cc,NCells,TR_local, tr, length(TRs{cc}));
        [BSpecMean_local, FitInfo_SpecMean]=lassoglm([XAmp XSpecMean],Y_local,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',1,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevSpecMean(tr),BestModSpecMean] = min(FitInfo_SpecMean.Deviance);
        BSpecMean(tr,2:end) = BSpecMean_local(:,BestModSpecMean);
        BSpecMean(tr,1) = FitInfo_SpecMean.Intercept(BestModSpecMean);
        
        % pitch saliency predicting Y
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d): Saliency Model\n', cc,NCells,TR_local, tr, length(TRs{cc}));
        [BSal_local, FitInfo_Sal]=lassoglm([XAmp XSaliency],Y_local,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',1,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevSal(tr),BestModSal] = min(FitInfo_Sal.Deviance);
        BSal(tr,2:end) = BSal_local(:,BestModSal);
        BSal(tr,1) = FitInfo_Sal.Intercept(BestModSal);
        
        
        %         % spectral median predicting Y
        %         [BSpecMed, FitInfo_SpecMed]=lassoglm(XSpecMed,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
        %         % find the model with the minimum of deviance (best lambda)
        %         [BestDevSpecMed(tr,dd),BestModSpecMed] = min(FitInfo_SpecMean.Deviance);
        
        % amplitude predicting Y, is a null model for the former 2 models
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d): Amplitude Model\n', cc,NCells,TR_local, tr, length(TRs{cc}));
        [BAmp_local, FitInfo_Amp]=lassoglm(XAmp,Y_local,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',1,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevAmp(tr),BestModAmp] = min(FitInfo_Amp.Deviance);
        
        BAmp(tr,2:end) = BAmp_local(:,BestModAmp);
        BAmp(tr,1) = FitInfo_Amp.Intercept(BestModAmp);
        
        % null model
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d): Null Model\n', cc,NCells,TR_local, tr, length(TRs{cc}));
        MDL_null=fitglm(ones(size(XAmp,1),1),Y_local,'Distribution',ParamModel.DISTR,'Link',ParamModel.LINK,'Intercept',false);
        % find the model with the minimum of deviance (best lambda)
        %     [BestDevNull(tr),BestModNull] = min(FitInfo_Amp.Deviance);
        %
        %     BNull(tr,2:end) = Bnull_local(:,BestModNull);
        %     BNull(tr,1) = FitInfo_null.Intercept(BestModNull);
        
        BestDevNull(tr) = MDL_null.Deviance;
        BNull(tr) = MDL_null.Coefficients.Estimate;
        
        if DatFig
            TimeBinsX = -Delay : (Win-Delay);
            TimeBinsX = TimeBinsX(2:end);
            figure(21)
            ColorCode = get(groot, 'DefaultAxesColorOrder');
            plot(TimeBinsX,BSpecMean_local((length(TimeBinsX)+1):end,BestModSpecMean),'LineStyle','-', 'LineWidth',2,'DisplayName', 'SpectralMean Betas','Color',ColorCode(1,:))
            hold on
            plot(TimeBinsX,BSpecMean_local(1:length(TimeBinsX),BestModSpecMean),'LineStyle',':', 'LineWidth',2,'DisplayName', 'SpectralMean Amp Betas','Color',ColorCode(1,:))
            hold on
            plot(TimeBinsX,BSal_local((length(TimeBinsX)+1):end,BestModSal),'LineStyle','-', 'LineWidth',2,'DisplayName','Saliency Betas' , 'Color',ColorCode(2,:))
            hold on
            plot(TimeBinsX,BSal_local(1:length(TimeBinsX),BestModSal),'LineStyle',':', 'LineWidth',2,'DisplayName','Saliency Amp Betas','Color',ColorCode(2,:))
            hold on
            %         plot(TimeBinsX,BSpecMed(:,BestModSpecMed),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Spectral Median Dev = %.1f',BestDevSpecMed(tr,dd)))
            hold on
            plot(TimeBinsX,BAmp_local(:,BestModAmp),'LineStyle','-', 'LineWidth',2,'DisplayName','Amp Betas' , 'Color', ColorCode(3,:))
            
            
            legend('show')
            % XTick = get(gca, 'XTickLabel');
            % XTick = cellfun(@str2double, XTick) * Win;
            % set(gca,'XTickLabel',XTick)
            xlabel('Time (ms)')
            title(sprintf('Model coefficients of Poisson Ridge regression on Acoustic Features, Delay = %dms Neural Response Time Resolution=%dms', Delay, TR_local))
            YLim = get(gca,'YLim');
            text(0,YLim(2)*0.8, sprintf('Spectral Mean Dev = %.1f',BestDevSpecMean(tr)))
            text(0,YLim(2)*0.7,sprintf('Saliency Dev = %.1f',BestDevSal(tr)));
            text(0,YLim(2)*0.6,sprintf('Amp Dev = %.1f',BestDevAmp(tr)))
            
            hold off
            pause(1)
        end
        TicToc(tr) = toc(TimerLoopGLM);
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d) => done in %.2f minutes \n', cc,NCells,TR_local, tr, length(TR_local), TicToc(tr)/60);
    end
    if DatFig
        
        figure(22)
        ColorCode = get(groot, 'DefaultAxesColorOrder');
        subplot(1,2,1)
        plot(TRs{cc},BestDevSpecMean','-o', 'MarkerFaceColor',ColorCode(1,:), 'Color', ColorCode(1,:),'LineWidth',2)
        hold on
        plot(TRs{cc},BestDevSal,'-o', 'MarkerFaceColor',ColorCode(2,:),'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs{cc},BestDevAmp,'-o', 'MarkerFaceColor',ColorCode(3,:), 'Color', ColorCode(3,:),'LineWidth',2)
        hold on
        plot(TRs{cc},BestDevNull,'-o', 'MarkerFaceColor',ColorCode(4,:), 'Color',ColorCode(4,:),'LineWidth',2, 'LineStyle','--')
        legend({'Amp + SpecMean' 'Amp + Sal' 'Amp' 'Null'})
        legend('Autoupdate','off','Location','southoutside')
        ylabel('Deviance')
        xlabel('Time resolution in ms')
        YLim = get(gca,'YLim');
        
        subplot(1,2,2)
        plot(TRs{cc},BestDevSpecMean', '-o', 'MarkerFaceColor',ColorCode(1,:),'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs{cc},BestDevSal,'-o', 'MarkerFaceColor',ColorCode(2,:),'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs{cc},BestDevAmp,'-o', 'MarkerFaceColor',ColorCode(3,:), 'Color', ColorCode(3,:),'LineWidth',2)
        hold on
        plot(TRs{cc},BestDevNull,'-o', 'MarkerFaceColor',ColorCode(4,:), 'Color',ColorCode(4,:),'LineWidth',2, 'LineStyle','--')
        ylabel('Deviance')
        xlabel('Time resolution in ms')
        xlim([0 2])
        ylim([0.9*YLim(2) YLim(2)])
        hold off
        
        
        figure(23)
        subplot(1,2,1)
        plot(TRs{cc}, BestDevAmp-BestDevSpecMean ,'-o', 'MarkerFaceColor',ColorCode(1,:), 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs{cc},BestDevAmp - BestDevSal,'-o', 'MarkerFaceColor',ColorCode(2,:),'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs{cc},BestDevNull -BestDevAmp,'-o', 'MarkerFaceColor',ColorCode(3,:), 'Color', ColorCode(3,:),'LineWidth',2)
        legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
        legend('Autoupdate','off','Location','southoutside')
        ylabel('Difference of Deviance')
        xlabel('Time resolution in ms')
        
        subplot(1,2,2)
        plot(TRs{cc}, (BestDevAmp-BestDevSpecMean)./BestDevAmp ,'-o', 'MarkerFaceColor',ColorCode(1,:), 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs{cc},(BestDevAmp - BestDevSal)./BestDevAmp,'-o', 'MarkerFaceColor',ColorCode(2,:),'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs{cc},(BestDevNull -BestDevAmp)./BestDevAmp,'-o', 'MarkerFaceColor',ColorCode(3,:), 'Color', ColorCode(3,:),'LineWidth',2)
        legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
        legend('Autoupdate','off','Location','southoutside')
        ylabel('R2')
        xlabel('Time resolution in ms')
        
        % Plotting betas
        figure(24)
        subplot(2,3,1:3)
        bar(TRs{cc},exp([BNull BAmp(:,1) BSpecMean(:,1) BSal(:,1)])*1000)
        xlabel('Time resolution in ms')
        ylabel('exp(intercept) Hz')
        legend({'Null' 'Amp' 'Amp + SpecMean' 'Amp + Sal'},'Location','southoutside','NumColumns',4)
        subplot(2,3,4)
        imagesc(BAmp(:,2:end))
        set(gca,'YTick',1:length(TRs{cc}),'YTickLabel',TRs{cc})
        ylabel('Time resolution in ms')
        xlabel('Time alligned to predicted Y (ms)')
        set(gca,'XTick',1:50:Win,'XTickLabel',-Delay:50:(Win-Delay-1))
        title('Amplitude Model')
        colorbar()
        subplot(2,3,5)
        imagesc(BSpecMean(:,2:end))
        set(gca,'YTick',1:length(TRs{cc}),'YTickLabel',TRs{cc})
        ylabel('Time resolution in ms')
        xlabel('Time alligned to predicted Y (ms)')
        set(gca,'XTick',1:100:2*Win,'XTickLabel',[-Delay:100:(Win-Delay-1) -Delay:100:(Win-Delay-1)])
        title('Amplitude + SpecMean')
        colorbar()
        subplot(2,3,6)
        imagesc(BSal(:,2:end))
        set(gca,'YTick',1:length(TRs{cc}),'YTickLabel',TRs{cc})
        ylabel('Time resolution in ms')
        xlabel('Time alligned to predicted Y (ms)')
        set(gca,'XTick',1:100:2*Win,'XTickLabel',[-Delay:100:(Win-Delay-1) -Delay:100:(Win-Delay-1)])
        title('Amplitude + Saliency')
        colorbar()
    end
    
    
    MotorModels{cc}.TRs = TRs{cc};
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
    
end
%
save(fullfile(Path,'MotorModelsGLM'),'MotorModels','CellsPath', 'GoodInfo');

%% Plot issues with no convergence and results of Poisson GLM
% load(fullfile(Path,'MotorModelsAllCells'),'MotorModels','CellsPath');
TicToc = nan(NCells*2,1);
MeanY = nan(NCells*2,1);
TRs = nan(NCells*2,1);
BestDevAmp = nan(NCells*2,1);
BestDevSpecMean = nan(NCells*2,1);
BestDevSal = nan(NCells*2,1);
BestDevNull = nan(NCells*2,1);

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
scatter(MeanY,TicToc/(60), 'MarkerFaceColor','k')
xlabel('MeanRate')
ylabel('Model Time consumption (min)')
subplot(1,2,2)
scatter(TRs,TicToc/(60), 'MarkerFaceColor','k')
xlabel('Models Time Resolution (ms)')
ylabel('Model Time consumption (min)')


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

% Now plot only results of the GLM for the time resolution dictated by
% width of time coherency peak
DevianceDifference = nan(NCells,3); % this matrix will contain the contribution of each parameter to the performance of the models predicting Y
R2_Models = nan(NCells,3); % this matrix will contain the contribution of each parameter to the performance of the models predicting Y

for cc=1:NCells
    if ~isempty(MotorModels{cc})
        DevianceDifference(cc,1) = MotorModels{cc}.DevNull(1) - MotorModels{cc}.DevAmp(1);
        DevianceDifference(cc,2) = MotorModels{cc}.DevAmp(1) - MotorModels{cc}.DevAmpSpecMean(1);
        DevianceDifference(cc,3) = MotorModels{cc}.DevAmp(1) - MotorModels{cc}.DevAmpSal(1);
        R2_Models(cc,1) = (MotorModels{cc}.DevNull(1) - MotorModels{cc}.DevAmp(1))/MotorModels{cc}.DevNull(1);
        R2_Models(cc,2) = (MotorModels{cc}.DevAmp(1) - MotorModels{cc}.DevAmpSpecMean(1))/MotorModels{cc}.DevAmp(1);
        R2_Models(cc,3) = (MotorModels{cc}.DevAmp(1) - MotorModels{cc}.DevAmpSal(1))/MotorModels{cc}.DevAmp(1);
    end
end

figure(30)
clf
subplot(1,2,1)
imagesc(DevianceDifference)
colorbar()
set(gca,'XTick',1:3, 'XTickLabel', {'Amp' 'SpecMean' 'Sal'})
ylabel('Cell #')
xlabel('Contribution to - Deviance')

subplot(1,2,2)
imagesc(R2_Models)
colorbar()
set(gca,'XTick',1:3, 'XTickLabel', {'Amp' 'SpecMean' 'Sal'})
%caxis([0 1])
colorbar()
ylabel('Cell #')
xlabel('Contribution to - Deviance in R2')

% Plot first cell example
figure(22)
clf
cc=1;
ColorCode = get(groot, 'DefaultAxesColorOrder');
X = categorical({ 'Null' 'Amp' 'Amp + SpecMean' 'Amp + Sal'});
X = reordercats(X,{ 'Null' 'Amp' 'Amp + SpecMean' 'Amp + Sal'});
bar(X,[MotorModels{cc}.DevNull(1) MotorModels{cc}.DevAmp(1) MotorModels{cc}.DevAmpSpecMean(1) MotorModels{cc}.DevAmpSal(1)],'k')
ylabel('Deviance = sum of errors')
title('Sum of errors (SSE for LM)')

figure(23)
clf
cc=1;
X = categorical({ 'Null - Amp' '(Amp + SpecMean) - Amp' '(Amp + Sal) - Amp'});
X = reordercats(X,{ 'Null - Amp' '(Amp + SpecMean) - Amp' '(Amp + Sal) - Amp'});
bar(X,[MotorModels{cc}.DevNull(1)-MotorModels{cc}.DevAmp(1) MotorModels{cc}.DevAmp(1)-MotorModels{cc}.DevAmpSpecMean(1) MotorModels{cc}.DevAmp(1)-MotorModels{cc}.DevAmpSal(1)],'k')
ylabel('Reduction in Deviance')
title('Absolute contribution of acoustic parameters to error reduction')

figure(24)
clf
cc=1;
ColorCode = get(groot, 'DefaultAxesColorOrder');
X = categorical({ 'Null - Amp' '(Amp + SpecMean) - Amp' '(Amp + Sal) - Amp'});
X = reordercats(X,{ 'Null - Amp' '(Amp + SpecMean) - Amp' '(Amp + Sal) - Amp'});
bar(X,[(MotorModels{cc}.DevNull(1)-MotorModels{cc}.DevAmp(1))/MotorModels{cc}.DevNull(1) (MotorModels{cc}.DevAmp(1)-MotorModels{cc}.DevAmpSpecMean(1))/MotorModels{cc}.DevAmp(1) (MotorModels{cc}.DevAmp(1)-MotorModels{cc}.DevAmpSal(1))/MotorModels{cc}.DevAmp(1)],'k')
ylabel('R2')
title('Relative contribution of acoustic parameters to error reduction')









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

function [XPerStim, XPerStimt] = get_x_4coherence(BioSound, Duration, Delay,TR,DefaultVal,Feature, Overlap)
DebugFig=0;
if nargin<7
    Overlap = 0;
end
if length(Delay)==1
    Delay = [Delay Delay];
end
% calculate the cell array of acoustic data for each cell. For each
% voc, XPerStim is a vector where each column correspond to the values of
% the acoustic feature in the window Win starting Delay ms before time onset
% of the vocalization.
%  Fs = round(1/((TR-Overlap).*10^-3));
% Vectors of the acoustic features
XPerStim = cell(1,length(Duration));
XPerStimt = cell(1,length(Duration));
if ischar(DefaultVal) && strcmp(DefaultVal, 'mean')
    DefaultValue = 1;
else
    DefaultValue = DefaultVal;
end

for stim = 1:length(Duration)
    if abs(round(length(BioSound{stim}.sound)/BioSound{stim}.samprate*10^3) - round(Duration(stim)))>1
        keyboard
    end
    % Get ready an output vector for the stim acoustic features that was sampled at 1000Hz
    XPerStim_temp = DefaultValue.*ones(1,round(1000 .* (Delay(1) + Duration(stim) + Delay(2)).*10^-3));
    FeatureVal = BioSound{stim}.(sprintf('%s',Feature));
    XPerStim_temp(Delay(1)+(1:length(FeatureVal)))=FeatureVal;
    if ~ischar(DefaultVal)
        XPerStim_temp(isnan(XPerStim_temp)) = DefaultValue;
    else
        XPerStim_temp(isnan(XPerStim_temp)) = nanmean(XPerStim_temp);
    end
    % Warning, resample gives very weird results here!!!
    %     [XPerStim{stim}] = resample(XPerStim_temp, Fs, 1000);
    %     [XPerStim{stim},XPerStimt{stim}] = resample(XPerStim_temp, (1:length(XPerStim_temp)),Fs*10^-3, 'spline');
    
    % I'm doing my own resampling
    % Time slots
    TimeBinsXOnsetInd = 1 :(TR-Overlap): (Delay(2) + Delay(1) + Duration(stim)); % These are slightly different than in get_Y_4GLM, because the times slot are used as indices in the vector and not as actuel time values!
    TimeBinsXOffsetInd = TimeBinsXOnsetInd + TR -1;
    TimeBinsXOnsetInd = TimeBinsXOnsetInd(TimeBinsXOffsetInd<=(Delay(2) + Delay(1) + Duration(stim))); % Only keep windows that are within the call
    TimeBinsXOffsetInd = TimeBinsXOffsetInd(TimeBinsXOffsetInd<=(Delay(2) + Delay(1) + Duration(stim))); % Only keep windows that are within the call
    XPerStim{stim} = nan(size(TimeBinsXOnsetInd));
    
    
    TimeBinsXOnset = -Delay(1) :(TR-Overlap): (Delay(2) + Duration(stim));
    TimeBinsXOffset = TimeBinsXOnset + TR;
    TimeBinsXOnset = TimeBinsXOnset(TimeBinsXOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    TimeBinsXOffset = TimeBinsXOffset(TimeBinsXOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    XPerStimt{stim} = TimeBinsXOnset + (TimeBinsXOffset - TimeBinsXOnset)/2;
    
    for tt=1:length(TimeBinsXOnsetInd)
        XPerStim{stim}(tt) = mean(XPerStim_temp(TimeBinsXOnsetInd(tt): TimeBinsXOffsetInd(tt)));
    end
    if DebugFig
        figure(200)
        clf
        plot((-Delay(1)+0.5):(Duration(stim)+Delay(2)),XPerStim_temp, 'LineWidth',2)
        xlabel('Time ms')
        ylabel(sprintf('%s',Feature))
        title(sprintf('Stim %d/%d',stim,length(Duration)));
        %         RemainTime = length(XPerStim_temp) - (length(XPerStim{stim})-1)*(1/Fs*10^3);
        %         XPerStimt{stim} = RemainTime/2+(1/Fs*10^3)*(0:(length(XPerStim{stim})-1));
        hold on
        plot(XPerStimt{stim},XPerStim{stim}, 'LineWidth',2)
        legend({'original' 'resampled'})
        pause(1)
    end
end


end


function [XPerStim,XPerStimt] = get_x_4GLM(BioSound, Duration, Delay,TR,DefaultVal,Feature, Overlap)
if nargin<7
    Overlap = 0;
end
if length(Delay)==1
    Delay = [Delay Delay];
end
% calculate the cell array of acoustic data for each cell. For each
% voc, XPerStim{stim} is a row vector where each column correspond to the value of
% the acoustic feature in the window TR starting Delay(1) ms before time onset
% of the vocalization.

% Vectors of the acoustic features
XPerStim = cell(1,length(Duration));
XPerStimt = cell(1,length(Duration));
if ischar(DefaultVal) && strcmp(DefaultVal, 'mean')
    DefaultValue = 1;
else
    DefaultValue = DefaultVal;
end

for stim = 1:length(Duration)
    if round(length(BioSound{stim}.sound)/BioSound{stim}.samprate*10^3) ~= Duration(stim)
        keyboard
    end
    
    % Time slots
    TimeBinsXOnsetInd = -Delay(1)+1 :(TR-Overlap): (Delay(2) + Duration(stim)); % These are slightly different than in get_Y_4GLM, because the times slot are used as indices in the vector and not as actuel time values!
    TimeBinsXOffsetInd = TimeBinsXOnsetInd + TR -1;
    TimeBinsXOnsetInd = TimeBinsXOnsetInd(TimeBinsXOffsetInd<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    TimeBinsXOffsetInd = TimeBinsXOffsetInd(TimeBinsXOffsetInd<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    XPerStim{stim} = nan(size(TimeBinsXOnsetInd));
    
    
    TimeBinsXOnset = -Delay(1) :(TR-Overlap): (Delay(2) + Duration(stim));
    TimeBinsXOffset = TimeBinsXOnset + TR;
    TimeBinsXOnset = TimeBinsXOnset(TimeBinsXOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    TimeBinsXOffset = TimeBinsXOffset(TimeBinsXOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    XPerStimt{stim} = TimeBinsXOnset + (TimeBinsXOffset - TimeBinsXOnset)/2;
    
    % Get ready an output vector for the stim acoustic features that was sampled at 1000Hz
    XPerStim_temp = DefaultValue.*ones(1,round(1000 .* (Delay(1) + Duration(stim) + Delay(2)).*10^-3));
    FeatureVal = BioSound{stim}.(sprintf('%s',Feature));
    XPerStim_temp(Delay(1)+(1:length(FeatureVal)))=FeatureVal;
    if ~ischar(DefaultVal)
        XPerStim_temp(isnan(XPerStim_temp)) = DefaultValue;
    else
        XPerStim_temp(isnan(XPerStim_temp)) = nanmean(XPerStim_temp);
    end
    for tt=1:length(TimeBinsXOnsetInd)
        XPerStim{stim}(tt) = mean(XPerStim_temp((TimeBinsXOnsetInd(tt)+Delay(1)): (TimeBinsXOffsetInd(tt)+Delay(1))));
        XPerStimt{stim}(tt) = mean(XPerStim_temp((TimeBinsXOnsetInd(tt)+Delay(1)): (TimeBinsXOffsetInd(tt)+Delay(1))));
    end
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

function [YPerStim, YPerStimt] = get_y_4Coherence(SAT, Duration,Delay,TR,Overlap)
DebugFig = 0;
if nargin<5
    Overlap = 0;
end
if length(Delay)==1
    Delay = [Delay Delay];
end
% Calculate the time varying rate applying a gaussian window TR on the
% spike pattern. The spike pattern considered starts -Delay ms
% before the onset of the vocalization and stops Delay ms after the
% offset of the vocalization
YPerStim = cell(1,length(Duration));
YPerStimt = cell(1,length(Duration));
% Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
nStd =(max(Duration) + Delay(1) + Delay(2))/10; % before set as 4
Tau = (TR/2);
T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
Expwav = Expwav./sum(Expwav);
% Frequency at which the neural data should be sampled
Fs = round(1/((TR-Overlap).*10^-3));
% Loop through the stimuli and fill in the matrix
for stim=1:length(Duration)
    % Time slots for the neural response
    TimeBinsY = -Delay(1) : (Delay(2) + Duration(stim));
    SpikePattern = zeros(1,length(TimeBinsY)-1);
    for isp = 1:length(SAT{stim})
        SpikeInd = round(SAT{stim}(isp));
        if (SpikeInd>=-Delay(1)) && (SpikeInd<(Delay(2) + Duration(stim)))
            SpikePattern(SpikeInd + Delay(1) +1) = SpikePattern(SpikeInd + Delay(1) +1) +1;
        end
    end
    
    % Convolve with Gaussian to obtain our smooth time varying spike train
    % and resample if necessary
    if Fs == 1000
        YPerStim{stim} = conv(SpikePattern, Expwav,'same');
        YPerStimt{stim} = TimeBinsY;
    else
        YPerStim_local = conv(SpikePattern, Expwav,'same');
        if sum(YPerStim_local)>0
            YPerStim_local = YPerStim_local/sum(YPerStim_local)*sum(SpikePattern); % Make sure we keep the right number of sipkes after convolution!
        end
        
        % resampling function is really doing weird things at edges...
        % doing my own resampling
        TimeBinsYOnsetInd = 1 :(TR-Overlap): (Delay(2) + Delay(1) + Duration(stim)); % These are slightly different than in get_Y_4GLM, because the times slot are used as indices in the vector and not as actuel time values!
        TimeBinsYOffsetInd = TimeBinsYOnsetInd + TR -1;
        TimeBinsYOnsetInd = TimeBinsYOnsetInd(TimeBinsYOffsetInd<=(Delay(2) + Delay(1) + Duration(stim))); % Only keep windows that are within the call
        TimeBinsYOffsetInd = TimeBinsYOffsetInd(TimeBinsYOffsetInd<=(Delay(2) + Delay(1) + Duration(stim))); % Only keep windows that are within the call
        YPerStim{stim} = nan(size(TimeBinsYOnsetInd));
        
        
        TimeBinsYOnset = -Delay(1) :(TR-Overlap): (Delay(2) + Duration(stim));
        TimeBinsYOffset = TimeBinsYOnset + TR;
        TimeBinsYOnset = TimeBinsYOnset(TimeBinsYOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
        TimeBinsYOffset = TimeBinsYOffset(TimeBinsYOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
        YPerStimt{stim} = TimeBinsYOnset + (TimeBinsYOffset - TimeBinsYOnset)/2;
        
        for tt=1:length(TimeBinsYOnsetInd)
            YPerStim{stim}(tt) = mean(YPerStim_local(TimeBinsYOnsetInd(tt): TimeBinsYOffsetInd(tt)));
        end
        if DebugFig
            figure(200)
            clf
            plot((-Delay(1)+0.5):(Duration(stim)+Delay(2)),YPerStim_local, 'LineWidth',2)
            xlabel('Time ms')
            ylabel('Spike Rate mHz (/ms)')
            title(sprintf('Stim %d/%d',stim,length(Duration)));
            %         RemainTime = length(XPerStim_temp) - (length(XPerStim{stim})-1)*(1/Fs*10^3);
            %         XPerStimt{stim} = RemainTime/2+(1/Fs*10^3)*(0:(length(XPerStim{stim})-1));
            hold on
            plot(YPerStimt{stim},YPerStim{stim}, 'LineWidth',2)
            legend({'original' 'resampled'}, 'AutoUpdate','off')
            hold on
            SpikeTimes = TimeBinsY(logical(SpikePattern));
            for ss = 1:length(SpikeTimes)
                V=vline(SpikeTimes(ss), 'k-');
                V.LineWidth = 2;
                hold on
            end
            pause(1)
        end
        
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
    
    % Make sure that the output mean(Y) = input mean(Y)
    %     if abs(round(sum(YPerStim{stim})*TR) - sum(SpikePattern))>TR/5
    %         warning('discrepancy in spike rate calculations larger than TR/5= %d?', TR/5)
    %         keyboard
    %     end
    
    %     if any(YPerStim{stim}<0)
    %         keyboard
    %     end
    
    
    %         % Time slots for the neural response
    %         TimeBinsY = -(Delay) : TR: (Delay + Duration(stim));
    %         YPerStim{stim} = nan(1,length(TimeBinsY)-1);
    %         for tt=1:(length(TimeBinsY)-1)
    %             % Find the number of spikes
    %             YPerStim{stim}(tt) = sum( (SAT{stim}>=TimeBinsY(tt)) .* (SAT{stim}<TimeBinsY(tt+1)));
    %         end
end

end







function [YPerStim] = get_y_4GLM(SAT, Duration,Delay,TR,Overlap)
if nargin<5
    Overlap = 0;
end
if length(Delay)==1
    Delay = [Delay Delay];
end
% Bin the spike arrival times in bins of TR ms with Overlap ms of overlap
% starting -Delay(1)ms before sound onset and stopping + Delay(2) ms after
% sound offset
YPerStim = cell(1,length(Duration));

% Loop through the stimuli and fill in the matrix
for stim=1:length(Duration)
    % Time slots for the neural response
    TimeBinsYOnset = -Delay(1) :(TR-Overlap): (Delay(2) + Duration(stim));
    TimeBinsYOffset = TimeBinsYOnset + TR;
    TimeBinsYOnset = TimeBinsYOnset(TimeBinsYOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    TimeBinsYOffset = TimeBinsYOffset(TimeBinsYOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    YPerStim{stim} = zeros(1,length(TimeBinsYOnset));
    % find the number of spikes for each window
    for tt = 1:length(TimeBinsYOnset)
        YPerStim{stim}(tt) = sum((SAT{stim}>=TimeBinsYOnset(tt)) .* (SAT{stim}<TimeBinsYOffset(tt)));
    end
    if any(YPerStim{stim}<0)
        keyboard
    end
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


function plotxyfeaturescoherence(BioSound,YPerStim,YPerStimt, XPerStim,XPerStimt,TR,Delay,Duration,F_high,FeatureName)
% This function plots for each stimulus the spectrogram with the
% corresponding acoustic feature and spike rate
if nargin<6
    FeatureName = 'Acoustic Feature';
end
if length(Delay)==1
    Delay = [Delay Delay];
end
DBNOISE =60;
f_low = 0;
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
    plot(XPerStimt{stim}, XPerStim{stim}, 'Color',ColorCode(2,:),'LineStyle','-', 'Marker','*', 'LineWidth',2)
    ylabel(sprintf('%s', FeatureName))
    if strcmp(FeatureName,'Spectral Mean')
        ylim(v_axis(3:4))
    elseif strcmp(FeatureName,'Saliency')
        YLim = get(gca, 'YLim');
        ylim([0 YLim(2)])
    elseif strcmp(FeatureName,'Amplitude')
        YLim = get(gca, 'YLim');
        ylim([0 YLim(2)])
    end
    
    xlim([-Delay(1) Duration(stim) + Delay(2)])
    hold off
    
    subplot(2,1,2)
    yyaxis left
    if length(YPerStim{stim})<4
        bar(YPerStimt{stim},YPerStim{stim}*10^3, 'BarWidth',0.2)
    else
        bar(YPerStimt{stim},YPerStim{stim}*10^3)
    end
    ylabel(sprintf('Spikes/s or Hz (%d ms Gauss win)', TR))
    YLIM = [0 500];
    ylim(YLIM)
    xlabel('Time (ms)')
    hold on
    line([0 Duration(stim)], [0.975 0.975].*YLIM(2), 'Color','k', 'LineWidth',12)
    hold on
    text(Duration(stim)/2.3,0.975.*YLIM(2),'Vocalization', 'Color', [1 1 1])
    hold on
    yyaxis right
    plot(XPerStimt{stim}, XPerStim{stim}, 'Color',ColorCode(2,:),'LineStyle','-', 'Marker','*', 'LineWidth',2)
    %     plot(TR/2 + (0:TR:((IndMax-1)*TR)), StimFeature(1:IndMax), 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    ylabel(sprintf('%s', FeatureName))
    if strcmp(FeatureName,'Spectral Mean')
        ylim(v_axis(3:4))
    elseif strcmp(FeatureName,'Saliency')
        YLim = get(gca, 'YLim');
        ylim([0 YLim(2)])
    elseif strcmp(FeatureName,'Amplitude')
        YLim = get(gca, 'YLim');
        ylim([0 YLim(2)])
    end
    xlim([-Delay(1) Duration(stim) + Delay(2)])
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


