addpath(genpath('/Users/elie/Documents/CODE/GitHub/SoundAnalysisBats'));
addpath(genpath('/Users/elie/Documents/CODE/GitHub/LoggerDataProcessing'));
addpath(genpath('/Users/elie/Documents/CODE/GitHub/LMC'));
addpath(genpath('/Users/elie/Documents/CODE/GitHub/GeneralCode'));
DatFig=0; %Set to 1 to see input data figures for each cell
OutFig = 1;%Set to 1 to see output data figures for each cell

%Filename = '59834_20190611_SSS_1-97.mat';
% Filename = '59834_20190610_SSS_1-130.mat';
%Path = '/Users/elie/Documents/LMCResults/';
HDPath = '/Volumes/JulieE8T/LMCResults/';
ServerPath = '/Volumes/server_home/users/JulieE/LMC/LMCResults';
% Path = '/Users/elie/Documents/ManipBats/LMC/ResultsFiles/';
% Path = '/Users/elie/Google Drive/BatmanData/';




%% List all cells
AllFiles = dir(fullfile(HDPath,'*8*.mat')); % Cells from both Cooper and Hodor
Files2run = zeros(length(AllFiles),1);
for ff=1:length(AllFiles)
    fprintf(1, 'Cell %d/%d\n', ff, length(AllFiles))
    Cell = load(fullfile(AllFiles(ff).folder, AllFiles(ff).name), 'What');
    if (length(strfind(AllFiles(ff).name, '_'))==3) && isfield(Cell, 'What') % Cells without the what field are files of another type or non-compiled cells probably because they are just noise
        Files2run(ff) = 1;
    end
    clear Cell
end
CellsPath = AllFiles(logical(Files2run));
NCells = length(CellsPath);

%% Check cell stability between 2 sessions
OpvsFr_Ks2stat = nan(NCells,1);
OpvsFr_Ks2statevenodd = nan(NCells,1);
OpvsFr_Ks2stat_boot = cell(NCells,1);
OpvsFr_Ks2pval_boot = nan(NCells,1);
OpvsFr_zs_KDE = nan(NCells,600); % contains the per cell normalized KDE (/min)
OpvsFr_meandiff_zs_KDE = nan(NCells,1);% mean difference zscored KDE between Operant and Free for each cell
Nbootstrap = 10000;
for cc=1:NCells
    
    %% load data
    CellPath = fullfile(CellsPath(cc).folder,CellsPath(cc).name);
    fprintf(1, '*** Cell %s %d/%d ***\n', CellPath, cc, NCells)
    Cell = load(CellPath,'QualitySSU', 'FreeBehavSession', 'OperantSession');
    % decide if cell is stable enough to compare free and operant sessions
    LogicalOp=logical(((Cell.QualitySSU.TimePoints(2:end)-60/2)>0) .* ((Cell.QualitySSU.TimePoints(2:end)-60/2)<diff(Cell.OperantSession)));
    LogicalFr=logical(((Cell.QualitySSU.TimePoints(2:end)-60/2)>(Cell.FreeBehavSession(1)- Cell.OperantSession(1))) .* ((Cell.QualitySSU.TimePoints(2:end)-60/2)<(Cell.FreeBehavSession(2)- Cell.OperantSession(1))));
    IndicesOp = find(LogicalOp);
    IndicesFr = find(LogicalFr);
    ZscoreKDE = (Cell.QualitySSU.KDE - nanmean(Cell.QualitySSU.KDE))/nanstd(Cell.QualitySSU.KDE);
    OpvsFr_zs_KDE(cc,(300-sum(LogicalOp)+1) : 300) = ZscoreKDE(LogicalOp);
    OpvsFr_zs_KDE(cc,300+(1:sum(LogicalFr))) = ZscoreKDE(LogicalFr);
    OpvsFr_meandiff_zs_KDE(cc) = mean(ZscoreKDE(LogicalOp)) - mean(ZscoreKDE(LogicalFr));
    if isempty(IndicesFr)
        fprintf(1, 'No Free session for this cell, no stability to calculate!\n')
        continue
    end
    AllIndices = [IndicesOp IndicesFr];
    EvenIndOp = IndicesOp(rem(IndicesOp,2)==0);
    OddIndOp = IndicesOp(rem(IndicesOp,2)==1);
    EvenIndFr = IndicesFr(rem(IndicesFr,2)==0);
    OddIndFr = IndicesFr(rem(IndicesFr,2)==1);
    %     [H,KDE_OpvsFr_p(cc),CI,STATS] = ttest2(Cell.QualitySSU.KDE(logical(IndicesOp)), Cell.QualitySSU.KDE(logical(IndicesFr)));
    if sum(LogicalFr)&& sum(LogicalOp)
        [~,~,OpvsFr_Ks2stat(cc)] = kstest2(Cell.QualitySSU.KDE(LogicalOp), Cell.QualitySSU.KDE(LogicalFr));
        [~,~,OpvsFr_Ks2statevenodd(cc)] = kstest2(Cell.QualitySSU.KDE([EvenIndOp EvenIndFr]), Cell.QualitySSU.KDE([OddIndOp OddIndFr]));
        OpvsFr_Ks2stat_boot_local = cell(1,Nbootstrap);
        parfor bb=1:Nbootstrap
            Shuffled = AllIndices(randperm(length(AllIndices)));
            [~,~,OpvsFr_Ks2stat_boot_local{bb}] = kstest2(Cell.QualitySSU.KDE(Shuffled(1:length(IndicesOp))), Cell.QualitySSU.KDE(Shuffled(length(IndicesOp)+ (1:length(IndicesFr)))));
        end
        OpvsFr_Ks2stat_boot{cc} = [OpvsFr_Ks2stat_boot_local{:}];
    end
    OpvsFr_Ks2pval_boot(cc) = sum(OpvsFr_Ks2stat_boot{cc}>=OpvsFr_Ks2stat(cc))/Nbootstrap;
    
%     figure(15);clf;subplot(1,2,1); histogram(Cell.QualitySSU.KDE(LogicalOp), 'Normalization', 'probability'); hold on; histogram(Cell.QualitySSU.KDE(LogicalFr), 'Normalization', 'probability')
%     xlabel('Spike Rate in Hz (bin=1min)')
%     ylabel('probability')
%     legend({'operant' 'Free'})
%     title('Spike rate distribution')
%     subplot(1,2,2)
%     histogram(OpvsFr_Ks2stat_boot{cc}, 'FaceColor', 'k'); hold on; VObs = vline(OpvsFr_Ks2stat(cc), 'g-'); VObs.LineWidth=2; hold on; VEO = vline(OpvsFr_Ks2statevenodd(cc), 'r-'); VEO.LineWidth=2;
%     title(sprintf('KS difference distribution p=%.2f',OpvsFr_Ks2pval_boot(cc)))
%     pause(1)
end
save(fullfile(HDPath, 'KDE_SessionStability.mat'), 'CellsPath', 'OpvsFr_Ks2stat', 'OpvsFr_Ks2statevenodd', 'OpvsFr_Ks2stat_boot','OpvsFr_Ks2pval_boot', 'OpvsFr_meandiff_zs_KDE', 'OpvsFr_zs_KDE')
save(fullfile(ServerPath, 'KDE_SessionStability.mat'), 'CellsPath', 'OpvsFr_Ks2stat', 'OpvsFr_Ks2statevenodd', 'OpvsFr_Ks2stat_boot','OpvsFr_Ks2pval_boot', 'OpvsFr_meandiff_zs_KDE', 'OpvsFr_zs_KDE')
figure();
subplot(1,3,1)
scatter(OpvsFr_Ks2stat-OpvsFr_Ks2statevenodd, OpvsFr_Ks2stat-cellfun(@mean, OpvsFr_Ks2stat_boot),30,[OpvsFr_Ks2pval_boot>0.01 zeros(size(OpvsFr_Ks2pval_boot)) zeros(size(OpvsFr_Ks2pval_boot))], 'filled'); 
xlabel('KS Obs - KS evenOdd'); ylabel('KS Obs - KS Bootstrap')
subplot(1,3,2)
LocalDat = OpvsFr_Ks2stat-cellfun(@mean, OpvsFr_Ks2stat_boot);
histogram(LocalDat(OpvsFr_Ks2pval_boot<=0.01), 'BinWidth',0.05); hold on; histogram(LocalDat(OpvsFr_Ks2pval_boot>0.01), 'BinWidth',0.05);xlabel('KS stat Obs - Bootstraped'); ylabel('# cells'); legend('Unstable','Stable')
subplot(1,3,3)
LocalDat = OpvsFr_Ks2stat-OpvsFr_Ks2statevenodd;
histogram(LocalDat(OpvsFr_Ks2pval_boot<=0.01), 'BinWidth',0.05); hold on; histogram(LocalDat(OpvsFr_Ks2pval_boot>0.01), 'BinWidth',0.05);xlabel('KS stat Obs - EvenOdd'); ylabel('# cells'); legend('Unstable','Stable')
suplabel(sprintf('Cell stability (N=%d): %d stable cells', NCells, sum(OpvsFr_Ks2pval_boot>0.01)))

FigRate = figure();
subplot(1,6,1:5)
[~,SortI] = sort(OpvsFr_meandiff_zs_KDE);
imagesc(OpvsFr_zs_KDE(SortI,1:500))
colormap("parula")
colorbar
Clim = caxis;
AmpC = diff(Clim);
caxis(Clim(1) + [0 AmpC*0.75])
xlabel('Time (min)')
ylabel('Cell #')
title('zscored-rate')
FigRate.CurrentAxes.XTick = [150 400];
FigRate.CurrentAxes.XTickLabel = {'Operant' 'Free'};
subplot(1,6,6)
scatter(OpvsFr_Ks2stat(SortI) - cellfun(@mean, OpvsFr_Ks2stat_boot(SortI)), 1:length(OpvsFr_Ks2stat), 20, [OpvsFr_Ks2pval_boot(SortI)>0.01 zeros(size(OpvsFr_Ks2pval_boot)) zeros(size(OpvsFr_Ks2pval_boot))], 'filled')
set(gca, 'YDir', 'reverse')
xlabel('Cell Instability (KS stat diff)')
%% Running through cells to find the optimal time resolution of the neural response for acoustic feature predicion from the neural response saving data per cell
% AllFiles = dir(fullfile(HDPath,'59834*.mat'));
AllFiles = dir(fullfile(HDPath,'11689*.mat'));
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
% BootstrapType=1;% 0:No Bootstrap; 1: bootstrap calculation based on Voc ID shuffling, keeping onsets; 2, shuffling in time of Y(neural response) within each voc keeping voc ID intact; 3, shuffling in time accross all neural vector used for coherence; 4:permutation of spikes accross all vocalizations
Average = 0; % 0: use feature variations for X; 1: use average feature variations accross vocalizations for X
NumVoc = 0;% 0: keep all vocalizations; any value: randomly choose NumVoc vocalizations to calculate Coherence
% Lags = -Delay:Delay;
% Freqs = (0:ceil(length(Lags)/2)).* (2*Nyquist/length(Lags)); % Lags is a uneven number so F(i) = i*2*Nyquist/length(Lags)


for cc=22:NCells
    
    CellTimer = tic();
    %% load data
    CellPath = fullfile(CellsPath(cc).folder,CellsPath(cc).name);
    fprintf(1, '*** Cell %s %d/%d %s ***\n', CellPath, cc, NCells, FeatureName)
    %     Cell = load(CellPath, 'What', 'ExpType', 'Who', 'BioSound', 'Duration','SpikesArrivalTimes_Behav','QualitySSU', 'VocOverlap', 'AudioQuality', 'DelayBefore','DelayAfter', 'AuditoryCoherenceFree', 'AuditoryCoherenceOperant', 'MotorCoherenceFree','MotorCoherenceOperant');
    Cell = load(CellPath, 'What', 'ExpType', 'Who', 'BioSound', 'Duration','SpikesArrivalTimes_Behav','QualitySSU', 'VocOverlap', 'AudioQuality', 'DelayBefore','DelayAfter');
    
    if ~Trill
        BootstrapType = 1; % bootstrap calculation based on Voc ID shuffling, keeping onsets
        AuditoryCoherenceFree = coherence4LMC(Cell,0,'Free', 'amp', Trill,NumVoc,BootstrapType, Average);
        MotorCoherenceFree = coherence4LMC(Cell, 1,'Free', 'amp', Trill,NumVoc,BootstrapType, Average);
        AuditoryCoherenceOperant = coherence4LMC(Cell, 0,'Operant', 'amp', Trill,NumVoc,BootstrapType, Average);
        MotorCoherenceOperant = coherence4LMC(Cell, 1,'Operant', 'amp', Trill,NumVoc,BootstrapType, Average);
        
        BootstrapType = 4; %permutation of spikes accross all vocalizations
        AuditoryCoherenceFree_local = coherence4LMC(Cell,0,'Free', 'amp', Trill,NumVoc,BootstrapType, Average);
        MotorCoherenceFree_local = coherence4LMC(Cell, 1,'Free', 'amp', Trill,NumVoc,BootstrapType, Average);
        AuditoryCoherenceOperant_local = coherence4LMC(Cell, 0,'Operant', 'amp', Trill,NumVoc,BootstrapType, Average);
        MotorCoherenceOperant_local = coherence4LMC(Cell, 1,'Operant', 'amp', Trill,NumVoc,BootstrapType, Average);
        
        
        if ~isfield(AuditoryCoherenceFree_local, 'Error')
            AuditoryCoherenceFree.Info_pRandSpikePerm = AuditoryCoherenceFree_local.Info_pRandSpikePerm;
            AuditoryCoherenceFree.BootstrapRandSpikePerm = AuditoryCoherenceFree_local.BootstrapRandSpikePerm;
        end
        
        if ~isfield(AuditoryCoherenceOperant_local, 'Error')
            AuditoryCoherenceOperant.Info_pRandSpikePerm = AuditoryCoherenceOperant_local.Info_pRandSpikePerm;
            AuditoryCoherenceOperant.BootstrapRandSpikePerm = AuditoryCoherenceOperant_local.BootstrapRandSpikePerm;
        end
        
        if ~isfield(MotorCoherenceFree_local, 'Error')
            MotorCoherenceFree.Info_pRandSpikePerm = MotorCoherenceFree_local.Info_pRandSpikePerm;
            MotorCoherenceFree.BootstrapRandSpikePerm = MotorCoherenceFree_local.BootstrapRandSpikePerm;
        end
        if ~isfield(MotorCoherenceOperant_local, 'Error')
            MotorCoherenceOperant.Info_pRandSpikePerm = MotorCoherenceOperant_local.Info_pRandSpikePerm;
            MotorCoherenceOperant.BootstrapRandSpikePerm = MotorCoherenceOperant_local.BootstrapRandSpikePerm;
        end
        fprintf(1, '\n*******************  saving results  ************************************\n')
        
%       save(CellPath, 'AuditoryCoherenceFree','MotorCoherenceFree','AuditoryCoherenceOperant', 'MotorCoherenceOperant','-append')
        
        warning('Error when saving to file, reading everything and rewriting over\n')
        Cell = load(CellPath);
        Cell.AuditoryCoherenceFree = AuditoryCoherenceFree;
        Cell.AuditoryCoherenceOperant = AuditoryCoherenceOperant;
        Cell.MotorCoherenceFree = MotorCoherenceFree;
        Cell.MotorCoherenceOperant = MotorCoherenceOperant;
        save(CellPath, '-struct','Cell')
        [~, FileName, Ext] = fileparts(CellPath);
        [S,~] = copyfile(CellPath, fullfile(ServerPath, [FileName Ext]));
        if ~S
            warning('File did not copy correcly to the server, advise!!')
            keyboard
        end
        
        
    elseif Trill==1 % Only calculating motor coherence for Trills because only Cooper is producing trills and never listening to them
        MotorCoherenceFreeTrill = coherence4LMC(Cell, 1,'Free', 'amp', Trill,NumVoc, BootstrapType);
        MotorCoherenceOperantTrill = coherence4LMC(Cell, 1,'Operant', 'amp', Trill,NumVoc, BootstrapType);
        if length(BootstrapType)==1 && (BootstrapType(1)==0) % only save to cell file if we are on the VocID shuffling for now
            save(CellPath, 'MotorCoherenceFreeTrill', 'MotorCoherenceOperantTrill','-append')
        end
        
        %         if ~isfield(MotorCoherenceFreeTrill, 'Error')
        %             MotorCoherenceFreeAll.CoherencyT_DelayAtzero(cc) = MotorCoherenceFreeTrill.CoherencyT_DelayAtzero;
        %             MotorCoherenceFreeAll.CoherencyT_WidthAtMaxPeak(cc) = MotorCoherenceFreeTrill.CoherencyT_WidthAtMaxPeak;
        %             MotorCoherenceFreeAll.MaxCoherence(cc) = MotorCoherenceFreeTrill.MaxCoherence1;
        %             MotorCoherenceFreeAll.MaxCoherenceF(cc) = MotorCoherenceFreeTrill.MaxCoherence2;
        %             MotorCoherenceFreeAll.CoherencePeaks{cc} = MotorCoherenceFreeTrill.CoherencePeaks;
        %             MotorCoherenceFreeAll.CoherencePeaksF{cc} = MotorCoherenceFreeTrill.CoherencePeaksF;
        %             MotorCoherenceFreeAll.FirstNonSigCoherenceFreq(cc) = MotorCoherenceFreeTrill.FirstNonSigCoherenceFreq;
        %             MotorCoherenceFreeAll.SecondCoherenceFreqCutOff(cc) = MotorCoherenceFreeTrill.SecondCoherenceFreqCutOff;
        %             MotorCoherenceFreeAll.Info(cc) = MotorCoherenceFreeTrill.Info;
        %             MotorCoherenceFreeAll.Info_low(cc) = MotorCoherenceFreeTrill.Info_low;
        %             MotorCoherenceFreeAll.Info_up(cc) = MotorCoherenceFreeTrill.Info_up;
        %             MotorCoherenceFreeAll.Info_p(cc) = MotorCoherenceFreeTrill.Info_p;
        %             MotorCoherenceFreeAll.Info_pTime(cc) = MotorCoherenceFreeTrill.Info_pTime;
        %         end
        %
        %
        %         if ~isfield(MotorCoherenceOperantTrill, 'Error')
        %             MotorCoherenceOperantAll.CoherencyT_DelayAtzero(cc) = MotorCoherenceOperantTrill.CoherencyT_DelayAtzero;
        %             MotorCoherenceOperantAll.CoherencyT_WidthAtMaxPeak(cc) = MotorCoherenceOperantTrill.CoherencyT_WidthAtMaxPeak;
        %             MotorCoherenceOperantAll.MaxCoherence(cc) = MotorCoherenceOperantTrill.MaxCoherence1;
        %             MotorCoherenceOperantAll.MaxCoherenceF(cc) = MotorCoherenceOperantTrill.MaxCoherence2;
        %             MotorCoherenceOperantAll.CoherencePeaks{cc} = MotorCoherenceOperantTrill.CoherencePeaks;
        %             MotorCoherenceOperantAll.CoherencePeaksF{cc} = MotorCoherenceOperantTrill.CoherencePeaksF;
        %             MotorCoherenceOperantAll.FirstNonSigCoherenceFreq(cc) = MotorCoherenceOperantTrill.FirstNonSigCoherenceFreq;
        %             MotorCoherenceOperantAll.SecondCoherenceFreqCutOff(cc) = MotorCoherenceOperantTrill.SecondCoherenceFreqCutOff;
        %             MotorCoherenceOperantAll.Info(cc) = MotorCoherenceOperantTrill.Info;
        %             MotorCoherenceOperantAll.Info_low(cc) = MotorCoherenceOperantTrill.Info_low;
        %             MotorCoherenceOperantAll.Info_up(cc) = MotorCoherenceOperantTrill.Info_up;
        %             MotorCoherenceOperantAll.Info_p(cc) = MotorCoherenceOperantTrill.Info_p;
        %             MotorCoherenceOperantAll.Info_pTime(cc) = MotorCoherenceOperantTrill.Info_pTime;
        %         end
    end
    
    
    fprintf(1, '******************** Done with Cell %s %d/%d in %ds ******************\n\n\n\n\n\n\n', CellPath,cc,NCells, toc(CellTimer));
    %             fprintf(1, 'Done with Cell %d/%d\n', cc,NCells)
    clear Cell
end


%% Collect values of coherence for plot purposes
AllFiles = dir(fullfile(HDPath,'*8*.mat')); % Cells from both Cooper and Hodor
Files2run = zeros(length(AllFiles),1);
FeatureName = 'amp';% {'amp' 'SpectralMean' 'sal'};
for ff=1:length(AllFiles)
    fprintf(1, 'Cell %d/%d\n', ff, length(AllFiles))
    Cell = load(fullfile(AllFiles(ff).folder, AllFiles(ff).name), 'What');
    if (length(strfind(AllFiles(ff).name, '_'))==3) && isfield(Cell, 'What') % Cells without the what field are files of another type or non-compiled cells probably because they are just noise
        Files2run(ff) = 1;
    end
    clear Cell
end
CellsPath = AllFiles(logical(Files2run));
NCells = length(CellsPath);
% initialize output variables
AuditoryCoherenceFreeAll = struct();
AuditoryCoherenceFreeAll.CoherencyT_DelayAtzero = nan(NCells,1);
AuditoryCoherenceFreeAll.CoherenceWeightedFreq = nan(NCells,1);
AuditoryCoherenceFreeAll.CumSumSigCoherence50Hz = nan(NCells,1);
AuditoryCoherenceFreeAll.MaxCoherence = nan(NCells,1);
AuditoryCoherenceFreeAll.MaxCoherenceF = nan(NCells,1);
AuditoryCoherenceFreeAll.CoherencePeaks = cell(NCells,1);
AuditoryCoherenceFreeAll.CoherencePeaksF = cell(NCells,1);
AuditoryCoherenceFreeAll.FirstNonSigCoherenceFreq = nan(NCells,1);
AuditoryCoherenceFreeAll.SecondCoherenceFreqCutOff = nan(NCells,1);
AuditoryCoherenceFreeAll.Info = nan(NCells,1);
AuditoryCoherenceFreeAll.Info_low = nan(NCells,1);
AuditoryCoherenceFreeAll.Info_up = nan(NCells,1);
AuditoryCoherenceFreeAll.Info_pRandSpikePerm = nan(NCells,1);
AuditoryCoherenceFreeAll.Info_pTime = nan(NCells,1);
AuditoryCoherenceFreeAll.Info_pRandVoc = nan(NCells,1);
AuditoryCoherenceFreeAll.Info_pFullTime = nan(NCells,1);
AuditoryCoherenceFreeAll.NStims = nan(NCells,1);
AuditoryCoherenceFreeAll.LengthX = nan(NCells,1);
AuditoryCoherenceFreeAll.AvAmpCoherencyT_DelayAtzero = nan(NCells,1);
AuditoryCoherenceFreeAll.AvAmpCoherencyT_WidthAtMaxPeak = nan(NCells,1);
AuditoryCoherenceFreeAll.AvAmpInfo = nan(NCells,1);
AuditoryCoherenceFreeAll.AvAmpInfo_low = nan(NCells,1);
AuditoryCoherenceFreeAll.AvAmpInfo_up = nan(NCells,1);
AuditoryCoherenceFreeAll.AvAmpInfo_pTime = nan(NCells,1);
AuditoryCoherenceFreeAll.CellsPath = CellsPath;
AuditoryCoherenceFreeAll.TR = TR;
AuditoryCoherenceFreeAll.Delay = Delay;
AuditoryCoherenceFreeAll.BatID = nan(NCells,1);

AuditoryCoherenceOperantAll = AuditoryCoherenceFreeAll;
MotorCoherenceFreeAll = AuditoryCoherenceFreeAll;
MotorCoherenceOperantAll = AuditoryCoherenceFreeAll;



for cc=1:NCells
    CellPath = fullfile(CellsPath(cc).folder,CellsPath(cc).name);
    fprintf(1, '*** Cell %s %d/%d %s ***\n', CellPath, cc, NCells, FeatureName)
    Cell = load(CellPath,  'AuditoryCoherenceFree', 'AuditoryCoherenceOperant', 'MotorCoherenceFree','MotorCoherenceOperant');
    if ~isfield(Cell.AuditoryCoherenceFree, 'Error')
        AuditoryCoherenceFreeAll.CoherencyT_DelayAtzero(cc) = Cell.AuditoryCoherenceFree.CoherencyT_DelayAtzero;
        AuditoryCoherenceFreeAll.CoherenceWeightedFreq(cc) = Cell.AuditoryCoherenceFree.WeightedSigCoherenceFreq;
        AuditoryCoherenceFreeAll.CumSumSigCoherence50Hz(cc) = Cell.AuditoryCoherenceFree.CumSumSigCoherence50Hz;
        AuditoryCoherenceFreeAll.MaxCoherence(cc) = Cell.AuditoryCoherenceFree.MaxCoherence1;
        AuditoryCoherenceFreeAll.MaxCoherenceF(cc) = Cell.AuditoryCoherenceFree.MaxCoherence2;
        AuditoryCoherenceFreeAll.CoherencePeaks{cc} = Cell.AuditoryCoherenceFree.CoherencePeaks;
        AuditoryCoherenceFreeAll.CoherencePeaksF{cc} = Cell.AuditoryCoherenceFree.CoherencePeaksF;
        AuditoryCoherenceFreeAll.FirstNonSigCoherenceFreq(cc) = Cell.AuditoryCoherenceFree.FirstNonSigCoherenceFreq;
        AuditoryCoherenceFreeAll.SecondCoherenceFreqCutOff(cc) = Cell.AuditoryCoherenceFree.SecondCoherenceFreqCutOff;
        AuditoryCoherenceFreeAll.Info(cc) = Cell.AuditoryCoherenceFree.Info;
        AuditoryCoherenceFreeAll.Info_low(cc) = Cell.AuditoryCoherenceFree.Info_low;
        AuditoryCoherenceFreeAll.Info_up(cc) = Cell.AuditoryCoherenceFree.Info_up;
        AuditoryCoherenceFreeAll.Info_pRandSpikePerm(cc) = Cell.AuditoryCoherenceFree.Info_pRandSpikePerm;
        AuditoryCoherenceFreeAll.Info_pRandVoc(cc) = Cell.AuditoryCoherenceFree.Info_pRandVoc;
        AuditoryCoherenceFreeAll.NStims(cc) = Cell.AuditoryCoherenceFree.NStims;
        AuditoryCoherenceFreeAll.LengthX(cc) = Cell.AuditoryCoherenceFree.LengthX;
        AuditoryCoherenceFreeAll.BatID(cc) = str2double(CellsPath(cc).name(1:5));
        
        %             AuditoryCoherenceFreeAll.Info_pTime(cc) = Cell.AuditoryCoherenceFree.Info_pTime;
        %             AuditoryCoherenceFreeAll.Info_pFullTime(cc) = Cell.AuditoryCoherenceFree.Info_pFullTime;
        %             AuditoryCoherenceFreeAll.AvAmpCoherencyT_DelayAtzero(cc) = Cell.AuditoryCoherenceFree_local.CoherencyT_DelayAtzero;
        %             AuditoryCoherenceFreeAll.AvAmpCoherencyT_WidthAtMaxPeak(cc) = Cell.AuditoryCoherenceFree_local.CoherencyT_WidthAtMaxPeak;
        %             AuditoryCoherenceFreeAll.AvAmpInfo(cc) = Cell.AuditoryCoherenceFree_local.Info;
        %             AuditoryCoherenceFreeAll.AvAmpInfo_low(cc) = Cell.AuditoryCoherenceFree_local.Info_low;
        %             AuditoryCoherenceFreeAll.AvAmpInfo_up(cc) = Cell.AuditoryCoherenceFree_local.Info_up;
        %             AuditoryCoherenceFreeAll.AvAmpInfo_pTime(cc) = sum((Cell.AuditoryCoherenceFree.BootstrapTime.Info - Cell.AuditoryCoherenceFree_local.Info)>= 0)/length(Cell.AuditoryCoherenceFree.BootstrapTime.Info);
        
    end
    if ~isfield(Cell.MotorCoherenceFree, 'Error')
        MotorCoherenceFreeAll.CoherencyT_DelayAtzero(cc) = Cell.MotorCoherenceFree.CoherencyT_DelayAtzero;
        MotorCoherenceFreeAll.CoherenceWeightedFreq(cc) = Cell.MotorCoherenceFree.WeightedSigCoherenceFreq;
        MotorCoherenceFreeAll.CumSumSigCoherence50Hz(cc) = Cell.MotorCoherenceFree.CumSumSigCoherence50Hz;
        MotorCoherenceFreeAll.MaxCoherence(cc) = Cell.MotorCoherenceFree.MaxCoherence1;
        MotorCoherenceFreeAll.MaxCoherenceF(cc) = Cell.MotorCoherenceFree.MaxCoherence2;
        MotorCoherenceFreeAll.CoherencePeaks{cc} = Cell.MotorCoherenceFree.CoherencePeaks;
        MotorCoherenceFreeAll.CoherencePeaksF{cc} = Cell.MotorCoherenceFree.CoherencePeaksF;
        MotorCoherenceFreeAll.FirstNonSigCoherenceFreq(cc) = Cell.MotorCoherenceFree.FirstNonSigCoherenceFreq;
        MotorCoherenceFreeAll.SecondCoherenceFreqCutOff(cc) = Cell.MotorCoherenceFree.SecondCoherenceFreqCutOff;
        MotorCoherenceFreeAll.Info(cc) = Cell.MotorCoherenceFree.Info;
        MotorCoherenceFreeAll.Info_low(cc) = Cell.MotorCoherenceFree.Info_low;
        MotorCoherenceFreeAll.Info_up(cc) = Cell.MotorCoherenceFree.Info_up;
        MotorCoherenceFreeAll.Info_pRandSpikePerm(cc) = Cell.MotorCoherenceFree.Info_pRandSpikePerm;
        MotorCoherenceFreeAll.Info_pRandVoc(cc) = Cell.MotorCoherenceFree.Info_pRandVoc;
        MotorCoherenceFreeAll.NStims(cc) = Cell.MotorCoherenceFree.NStims;
        MotorCoherenceFreeAll.LengthX(cc) = Cell.MotorCoherenceFree.LengthX;
        MotorCoherenceFreeAll.BatID(cc) = str2double(CellsPath(cc).name(1:5));
        
        %             MotorCoherenceFreeAll.Info_pTime(cc) = Cell.MotorCoherenceFree.Info_pTime;
        %             MotorCoherenceFreeAll.Info_pFullTime(cc) = Cell.MotorCoherenceFree.Info_pFullTime;
        %             MotorCoherenceFreeAll.AvAmpCoherencyT_DelayAtzero(cc) = Cell.MotorCoherenceFree_local.CoherencyT_DelayAtzero;
        %             MotorCoherenceFreeAll.AvAmpCoherencyT_WidthAtMaxPeak(cc) = Cell.MotorCoherenceFree_local.CoherencyT_WidthAtMaxPeak;
        %             MotorCoherenceFreeAll.AvAmpInfo(cc) = Cell.MotorCoherenceFree_local.Info;
        %             MotorCoherenceFreeAll.AvAmpInfo_low(cc) = Cell.MotorCoherenceFree_local.Info_low;
        %             MotorCoherenceFreeAll.AvAmpInfo_up(cc) = Cell.MotorCoherenceFree_local.Info_up;
        %             MotorCoherenceFreeAll.AvAmpInfo_pTime(cc) = sum((Cell.MotorCoherenceFree.BootstrapTime.Info - Cell.MotorCoherenceFree_local.Info)>= 0)/length(Cell.MotorCoherenceFree.BootstrapTime.Info);
    end
    
    if ~isfield(Cell.AuditoryCoherenceOperant, 'Error')
        AuditoryCoherenceOperantAll.CoherencyT_DelayAtzero(cc) = Cell.AuditoryCoherenceOperant.CoherencyT_DelayAtzero;
        AuditoryCoherenceOperantAll.CoherenceWeightedFreq(cc) = Cell.AuditoryCoherenceOperant.WeightedSigCoherenceFreq;
        AuditoryCoherenceOperantAll.CumSumSigCoherence50Hz(cc) = Cell.AuditoryCoherenceOperant.CumSumSigCoherence50Hz;
        AuditoryCoherenceOperantAll.MaxCoherence(cc) = Cell.AuditoryCoherenceOperant.MaxCoherence1;
        AuditoryCoherenceOperantAll.MaxCoherenceF(cc) = Cell.AuditoryCoherenceOperant.MaxCoherence2;
        AuditoryCoherenceOperantAll.CoherencePeaks{cc} = Cell.AuditoryCoherenceOperant.CoherencePeaks;
        AuditoryCoherenceOperantAll.CoherencePeaksF{cc} = Cell.AuditoryCoherenceOperant.CoherencePeaksF;
        AuditoryCoherenceOperantAll.FirstNonSigCoherenceFreq(cc) = Cell.AuditoryCoherenceOperant.FirstNonSigCoherenceFreq;
        AuditoryCoherenceOperantAll.SecondCoherenceFreqCutOff(cc) = Cell.AuditoryCoherenceOperant.SecondCoherenceFreqCutOff;
        AuditoryCoherenceOperantAll.Info(cc) = Cell.AuditoryCoherenceOperant.Info;
        AuditoryCoherenceOperantAll.Info_low(cc) = Cell.AuditoryCoherenceOperant.Info_low;
        AuditoryCoherenceOperantAll.Info_up(cc) = Cell.AuditoryCoherenceOperant.Info_up;
        AuditoryCoherenceOperantAll.Info_pRandSpikePerm(cc) = Cell.AuditoryCoherenceOperant.Info_pRandSpikePerm;
        AuditoryCoherenceOperantAll.Info_pRandVoc(cc) = Cell.AuditoryCoherenceOperant.Info_pRandVoc;
        AuditoryCoherenceOperantAll.NStims(cc) = Cell.AuditoryCoherenceOperant.NStims;
        AuditoryCoherenceOperantAll.LengthX(cc) = Cell.AuditoryCoherenceOperant.LengthX;
        AuditoryCoherenceOperantAll.BatID(cc) = str2double(CellsPath(cc).name(1:5));
    end
    
    if ~isfield(Cell.MotorCoherenceOperant, 'Error')
        MotorCoherenceOperantAll.CoherencyT_DelayAtzero(cc) = Cell.MotorCoherenceOperant.CoherencyT_DelayAtzero;
        MotorCoherenceOperantAll.CoherenceWeightedFreq(cc) = Cell.MotorCoherenceOperant.WeightedSigCoherenceFreq;
        MotorCoherenceOperantAll.CumSumSigCoherence50Hz(cc) = Cell.MotorCoherenceOperant.CumSumSigCoherence50Hz;
        MotorCoherenceOperantAll.MaxCoherence(cc) = Cell.MotorCoherenceOperant.MaxCoherence1;
        MotorCoherenceOperantAll.MaxCoherenceF(cc) = Cell.MotorCoherenceOperant.MaxCoherence2;
        MotorCoherenceOperantAll.CoherencePeaks{cc} = Cell.MotorCoherenceOperant.CoherencePeaks;
        MotorCoherenceOperantAll.CoherencePeaksF{cc} = Cell.MotorCoherenceOperant.CoherencePeaksF;
        MotorCoherenceOperantAll.FirstNonSigCoherenceFreq(cc) = Cell.MotorCoherenceOperant.FirstNonSigCoherenceFreq;
        MotorCoherenceOperantAll.SecondCoherenceFreqCutOff(cc) = Cell.MotorCoherenceOperant.SecondCoherenceFreqCutOff;
        MotorCoherenceOperantAll.Info(cc) = Cell.MotorCoherenceOperant.Info;
        MotorCoherenceOperantAll.Info_low(cc) = Cell.MotorCoherenceOperant.Info_low;
        MotorCoherenceOperantAll.Info_up(cc) = Cell.MotorCoherenceOperant.Info_up;
        MotorCoherenceOperantAll.Info_pRandVoc(cc) = Cell.MotorCoherenceOperant.Info_pRandVoc;
        MotorCoherenceOperantAll.Info_pRandSpikePerm(cc) = Cell.MotorCoherenceOperant.Info_pRandSpikePerm;
        MotorCoherenceOperantAll.NStims(cc) = Cell.MotorCoherenceOperant.NStims;
        MotorCoherenceOperantAll.LengthX(cc) = Cell.MotorCoherenceOperant.LengthX;
        MotorCoherenceOperantAll.BatID(cc) = str2double(CellsPath(cc).name(1:5));
        
        %             MotorCoherenceOperantAll.Info_pTime(cc) = Cell.MotorCoherenceOperant.Info_pTime;
        %             MotorCoherenceOperantAll.Info_pFullTime(cc) = Cell.MotorCoherenceOperant.Info_pFullTime;
        %             MotorCoherenceOperantAll.AvAmpCoherencyT_DelayAtzero(cc) = Cell.MotorCoherenceOperant_local.CoherencyT_DelayAtzero;
        %             MotorCoherenceOperantAll.AvAmpCoherencyT_WidthAtMaxPeak(cc) = Cell.MotorCoherenceOperant_local.CoherencyT_WidthAtMaxPeak;
        %             MotorCoherenceOperantAll.AvAmpInfo(cc) = Cell.MotorCoherenceOperant_local.Info;
        %             MotorCoherenceOperantAll.AvAmpInfo_low(cc) = Cell.MotorCoherenceOperant_local.Info_low;
        %             MotorCoherenceOperantAll.AvAmpInfo_up(cc) = Cell.MotorCoherenceOperant_local.Info_up;
        %             MotorCoherenceOperantAll.AvAmpInfo_pTime(cc) = sum((Cell.MotorCoherenceOperant.BootstrapTime.Info - Cell.MotorCoherenceOperant_local.Info)>= 0)/length(Cell.MotorCoherenceOperant.BootstrapTime.Info);
    end
    
end


% Order cells by
% decreasing values of info
if ~Trill
    [~,AuditoryCoherenceFreeAll.GoodInfo] = sort(AuditoryCoherenceFreeAll.Info,'descend');
    save(fullfile(HDPath,sprintf('AuditoryCoherence_%s_%s.mat', FeatureName, 'Free')),'-struct','AuditoryCoherenceFreeAll');
    save(fullfile(ServerPath,sprintf('AuditoryCoherence_%s_%s.mat', FeatureName, 'Free')),'-struct','AuditoryCoherenceFreeAll');
    
    [~,MotorCoherenceFreeAll.GoodInfo] = sort(MotorCoherenceFreeAll.Info,'descend');
    save(fullfile(HDPath,sprintf('MotorCoherence_%s_%s.mat', FeatureName, 'Free')),'-struct','MotorCoherenceFreeAll');
    save(fullfile(ServerPath,sprintf('MotorCoherence_%s_%s.mat', FeatureName, 'Free')),'-struct','MotorCoherenceFreeAll');
    
    [~,AuditoryCoherenceOperantAll.GoodInfo] = sort(AuditoryCoherenceOperantAll.Info,'descend');
    save(fullfile(HDPath,sprintf('AuditoryCoherence_%s_%s.mat', FeatureName, 'Operant')),'-struct','AuditoryCoherenceOperantAll');
    save(fullfile(ServerPath,sprintf('AuditoryCoherence_%s_%s.mat', FeatureName, 'Operant')),'-struct','AuditoryCoherenceOperantAll');
    
    [~,MotorCoherenceOperantAll.GoodInfo] = sort(MotorCoherenceOperantAll.Info,'descend');
    save(fullfile(HDPath,sprintf('MotorCoherence_%s_%s.mat', FeatureName, 'Operant')),'-struct','MotorCoherenceOperantAll');
    save(fullfile(ServerPath,sprintf('MotorCoherence_%s_%s.mat', FeatureName, 'Operant')),'-struct','MotorCoherenceOperantAll');
elseif Trill
    [~,MotorCoherenceFreeAll.GoodInfo] = sort(MotorCoherenceFreeAll.Info,'descend');
    save(fullfile(HDPath,sprintf('MotorCoherence_%s_%s_Trill.mat', FeatureName, 'Free')),'-struct','MotorCoherenceFreeAll');
    
    [~,MotorCoherenceOperantAll.GoodInfo] = sort(MotorCoherenceOperantAll.Info,'descend');
    save(fullfile(HDPath,sprintf('MotorCoherence_%s_%s_Trill.mat', FeatureName, 'Operant')),'-struct','MotorCoherenceOperantAll');
else
    warning('Code not written for these conditions!')
    keyboard
end

%% Plots results of coherence calculations for the population
FeatureName = {'amp' 'SpectralMean' 'sal'};
SessionChoice = 'Operant'; % Operant or Free
MS = 'Motor'; % Motor or Auditory
Trill=0;
% for fn = 1:length(FeatureName)
fn=1;
if ~Trill
    load(fullfile(HDPath,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, SessionChoice)))
elseif Trill==1
    load(fullfile(HDPath,sprintf('%sCoherence_%s_%s_Trill.mat', MS, FeatureName{fn}, SessionChoice)))
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
Colors = [(Info_pRandSpikePerm<0.01) zeros(length(Info),1) zeros(length(Info),1)];
scatter(Info(BatID==59834),MaxCoherence(BatID==59834,1),20, Colors(BatID==59834,:), 'filled')
hold on
scatter(Info(BatID==11689),MaxCoherence(BatID==11689,1),20, Colors(BatID==11689,:), 'filled', '^')
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
    suplabel(sprintf('%s Coherence %s %s session (N = %d + %d = %d)', MS, FeatureName2, SessionChoice, sum(~isnan(Info(BatID==59834))),sum(~isnan(Info(BatID==11689))),sum(~isnan(Info))), 't');
elseif Trill==1
    suplabel(sprintf('%s Coherence %s Trills only %s session (N = %d)', MS, FeatureName2, SessionChoice, sum(~isnan(Info))), 't');
end


% Just the histogram of information
F3=figure(3);
clf
histogram(Info,'BinWidth',0.1,'FaceColor', [0 0.447 0.741])
xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
ylabel('Number of Cells')
hold on
v=vline(nanmean(Info),'c--');
v.LineWidth = 2;
XLIM = get(gca, 'XLim');
set(gca, 'XLim', [XLIM(1) 5])
title(sprintf('%s Information (mean = %.2f +/- %.2f, N=%d)', MS, nanmean(Info), nanstd(Info), sum(~isnan(Info))))
box off
hold off

% Just the histogram of information with cell breaked down
GestureCell_logical = logical((Info_pRandSpikePerm<0.01).*(Info_pRandVoc<0.01));% Gesture cells
OnCell_logical = logical((Info_pRandSpikePerm<0.01).*(Info_pRandVoc>=0.01));% On cells
NSCell_logical = Info_pRandSpikePerm>=0.01;
SCell_logical = Info_pRandSpikePerm<0.01;
FBreak=figure(3);
clf
histogram(Info(OnCell_logical),'BinWidth',0.1,'FaceColor', [0, 0, 1]) % Vocalization On Cells; purple
hold on
histogram(Info(GestureCell_logical),'BinWidth',0.1,'FaceColor', [0.75, 0, 0.75]) % Gesture Cells only; red
hold on
histogram(Info(NSCell_logical),'BinWidth',0.1,'FaceColor', 'k') % Non Significant cells
legend({'Vocalization On units', 'Gesture units', 'NS units'})
xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
ylabel('Number of Cells')
hold on
v=vline(nanmean(Info(~NSCell_logical)),'c--');
v.LineWidth = 2;
XLIM = get(gca, 'XLim');
set(gca, 'XLim', [XLIM(1) 5])
title(sprintf('%s Information (mean = %.2f +/- %.2f, N=%d)', MS, nanmean(Info(SCell_logical)), nanstd(Info(SCell_logical)),sum(SCell_logical)))
box off
hold off

% Compile figure for delay and time resolution
figure(4)
clf
t = tiledlayout(1,1);
ax1 = axes(t);
histogram(ax1,CoherencyT_DelayAtzero(SCell_logical),'BinWidth',2*TR,'FaceColor', [0.85 0.325 0.098], 'orientation', 'horizontal')

%     hold(ax1, 'on')
%     histogram(CoherencyT_DelayAtzero(Info_p<0.01),'BinWidth',TR,'FaceColor','r')

xlabel('Number of cells')
ax1.XAxisLocation = 'top';
ax1.YAxisLocation = 'right';
ax1.Box = 'off';
ax1.Color = 'none';
text(2,-450, sprintf('Sig mean = %.2f +/- %.2f ms', nanmean(CoherencyT_DelayAtzero(SCell_logical)), nanstd(CoherencyT_DelayAtzero(SCell_logical))/(sum(~isnan(CoherencyT_DelayAtzero(SCell_logical))))^0.5))
%     hold on
%     text(29, -57, sprintf('Sig mean = %.2f +/- %.2f', nanmean(CoherencyT_DelayAtzero(Info_p<0.01)), nanstd(CoherencyT_DelayAtzero(Info_p<0.01))/(sum(~isnan(CoherencyT_DelayAtzero(Info_p<0.01))))^0.5))
hold(ax1,'on')
ax2 = axes(t);
scatter(ax2,Info(logical(OnCell_logical.*(BatID == 59834))),CoherencyT_DelayAtzero(logical(OnCell_logical.*(BatID == 59834))),25, [0, 0, 1], 'filled')
hold(ax2, 'on')
scatter(ax2,Info(logical(OnCell_logical.*(BatID == 11689))),CoherencyT_DelayAtzero(logical(OnCell_logical.*(BatID == 11689))),25, [0, 0, 1], 'filled','d')
hold(ax2, 'on')
scatter(ax2,Info(logical(GestureCell_logical.*(BatID == 59834))),CoherencyT_DelayAtzero(logical(GestureCell_logical.*(BatID == 59834))),25, [0.75, 0, 0.75], 'filled')
hold(ax2, 'on')
scatter(ax2,Info(logical(GestureCell_logical.*(BatID == 11689))),CoherencyT_DelayAtzero(logical(GestureCell_logical.*(BatID == 11689))),25, [0.75, 0, 0.75], 'filled','d')
%     hold(ax2, 'on')
%     scatter(ax2,Info(NSCell_logical),CoherencyT_DelayAtzero(NSCell_logical),25, 'k', 'filled')
legend({'Cooper Vocalization On units','Hodor Vocalization On units', 'Cooper Gesture units', 'Hodor Gesture units'})
%     legend({'Vocalization On units', 'Gesture units', 'NS units'})
ylabel('Optimal Delay (ms)')
ax2.Box = 'off';
ax2.Color = 'none';
xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
xlim([0 5])
hold on
h=hline(0, 'k:');
h.LineWidth = 2;
ylim(ax1, ylim(ax2));

OptimalTimeResolution = 1000./(2*CoherenceWeightedFreq);
figure(5)
clf
t = tiledlayout(1,1);
ax1 = axes(t);
histogram(ax1, OptimalTimeResolution(SCell_logical),'BinWidth',TR,'FaceColor', [0 0.5 0], 'orientation', 'horizontal')
%     ylabel('Time resolution (ms)')
xlabel('Number of cells')
ax1.XAxisLocation = 'top';
ax1.YAxisLocation = 'right';
ax1.Box = 'off';
ax1.Color = 'none';
ax1.XLim = xlim(ax1) +[0 1];
ax1.YLim = [0 45];
text(20, 3, sprintf('Sig mean = %.2f +/- %.2f', nanmean(OptimalTimeResolution(SCell_logical)), nanstd(OptimalTimeResolution(SCell_logical))/(sum(~isnan(OptimalTimeResolution(SCell_logical))))^0.5))
%     hold on
%     text(-57, 29, sprintf('Sig mean = %.2f +/- %.2f', nanmean(OptimalTimeResolution(Info_p<0.01)), nanstd(OptimalTimeResolution(Info_p<0.01))/(sum(~isnan(OptimalTimeResolution(Info_p<0.01))))^0.5))

hold(ax1,'on')
ax2 = axes(t);
scatter(ax2,Info(logical(OnCell_logical.*(BatID == 59834))),OptimalTimeResolution(logical(OnCell_logical.*(BatID == 59834))),25, [0, 0, 1], 'filled')
hold(ax2, 'on')
scatter(ax2,Info(logical(OnCell_logical.*(BatID == 11689))),OptimalTimeResolution(logical(OnCell_logical.*(BatID == 11689))),25, [0, 0, 1], 'filled','d')
hold(ax2, 'on')
scatter(ax2,Info(logical(GestureCell_logical.*(BatID == 59834))),OptimalTimeResolution(logical(GestureCell_logical.*(BatID == 59834))),25, [0.75, 0, 0.75], 'filled')
hold(ax2, 'on')
scatter(ax2,Info(logical(GestureCell_logical.*(BatID == 11689))),OptimalTimeResolution(logical(GestureCell_logical.*(BatID == 11689))),25, [0.75, 0, 0.75], 'filled','d')
%     hold(ax2, 'on')
%     scatter(ax2,Info(NSCell_logical),OptimalTimeResolution(NSCell_logical),25, 'k', 'filled')
%     legend({'Vocalization On units', 'Gesture units', 'NS units'})
legend({'Cooper Vocalization On units', 'Hodor Vocalization On units', 'Cooper Gesture units', 'Hodor Gesture units'})
ax2.Box = 'off';
ax2.Color = 'none';
xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
ylabel('Time resolution (ms)')
xlim([0 5])
ax2.YLim = ylim(ax1);




figure(6)
clf
subplot(2,5,1)
Colors = [(Info_pRandSpikePerm<0.01) zeros(length(Info),1) zeros(length(Info),1)];
scatter(Info(BatID==59834),CoherencyT_DelayAtzero(BatID==59834),20, Colors(BatID==59834,:), 'filled')
hold on
scatter(Info(BatID==11689),CoherencyT_DelayAtzero(BatID==11689),20, Colors(BatID==11689,:), 'filled', 'd')
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
scatter(Info(BatID==59834),OptimalTimeResolution(BatID==59834),20, Colors(BatID==59834,:), 'filled')
hold on
scatter(Info(BatID==11689),OptimalTimeResolution(BatID==11689),20, Colors(BatID==11689,:), 'filled', 'd')
xlabel(sprintf('Information on Coherence with sound %s (bits)',FeatureName2))
ylabel('Time resolution of Coherence (ms)')
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
scatter(Info(BatID==59834),FirstNonSigCoherenceFreq(BatID==59834),20, Colors(BatID==59834,:), 'filled')
hold on
scatter(Info(BatID==11689),FirstNonSigCoherenceFreq(BatID==11689),20, Colors(BatID==11689,:), 'filled','d')
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
scatter(Info(BatID==59834),SecondCoherenceFreqCutOff(BatID==59834),20, Colors(BatID==59834,:), 'filled')
hold on
scatter(Info(BatID==11689),SecondCoherenceFreqCutOff(BatID==11689),20, Colors(BatID==11689,:), 'filled', 'd')
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
scatter(Info(BatID==59834),MaxCoherenceF(BatID==59834),20, Colors(BatID==59834,:), 'filled')
hold on
scatter(Info(BatID==11689),MaxCoherenceF(BatID==11689),20, Colors(BatID==11689,:), 'filled')
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
histogram(CoherencyT_DelayAtzero(Info_pRandSpikePerm>=0.01),'BinWidth',TR/2,'FaceColor','k')
hold on
histogram(CoherencyT_DelayAtzero(Info_pRandSpikePerm<0.01),'BinWidth',TR/2,'FaceColor','r')
xlabel('Phase of coherency (ms)')
ylabel('Number of cells')
hold on
v=vline(0, 'r:');
v.LineWidth = 2;
hold on
text(-57, 27, sprintf('All mean = %.2f +/- %.2f', nanmean(CoherencyT_DelayAtzero), nanstd(CoherencyT_DelayAtzero)/(sum(~isnan(CoherencyT_DelayAtzero)))^0.5))
hold on
text(-57, 29, sprintf('Sig mean = %.2f +/- %.2f', nanmean(CoherencyT_DelayAtzero(Info_pRandSpikePerm<0.01)), nanstd(CoherencyT_DelayAtzero(Info_pRandSpikePerm<0.01))/(sum(~isnan(CoherencyT_DelayAtzero(Info_pRandSpikePerm<0.01))))^0.5))

subplot(2,5,7)
histogram(OptimalTimeResolution(Info_pRandSpikePerm>=0.01),'BinWidth',TR/2,'FaceColor','k')
hold on
histogram(OptimalTimeResolution(Info_pRandSpikePerm<0.01),'BinWidth',TR/2,'FaceColor','r')
xlabel('Time resolution of Coherency (ms)')
ylabel('Number of cells')
hold on
text(35, 32, sprintf('All mean = %.2f +/- %.2f', nanmean(OptimalTimeResolution), nanstd(OptimalTimeResolution)/(sum(~isnan(OptimalTimeResolution)))^0.5))
hold on
text(35, 34, sprintf('Sig mean = %.2f +/- %.2f', nanmean(OptimalTimeResolution(Info_pRandSpikePerm<0.01)), nanstd(OptimalTimeResolution(Info_pRandSpikePerm<0.01))/(sum(~isnan(OptimalTimeResolution(Info_pRandSpikePerm<0.01))))^0.5))



subplot(2,5,8)
histogram(FirstNonSigCoherenceFreq(Info_pRandSpikePerm>=0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
hold on
histogram(FirstNonSigCoherenceFreq(Info_pRandSpikePerm<0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','r')
xlabel('Max significant Frequency (Hz)')
ylabel('Number of Cells')

subplot(2,5,9)
histogram(SecondCoherenceFreqCutOff(Info_pRandSpikePerm>=0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
hold on
histogram(SecondCoherenceFreqCutOff(Info_pRandSpikePerm<0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','r')
xlabel('Max 2nd peaks significant Frequency (Hz)')
ylabel('Number of Cells')

subplot(2,5,10)
histogram(MaxCoherenceF(Info_pRandSpikePerm>=0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','k')
hold on
histogram(MaxCoherenceF(Info_pRandSpikePerm<0.01), 'BinWidth',ceil(nFFT/Nyquist),'FaceColor','r')
% Some points with very low values of info have high values of frequency of
% Max coherence, keep the plot focused on the majority of points
xlim([0 50])
xlabel('Frequency of Max Coherence (Hz)')
ylabel('Number of Cells')

if ~Trill
    suplabel(sprintf('%s Coherence %s session %s N=%d cells %d with spike shuffling test p<0.01', MS, SessionChoice, FeatureName2, sum(~isnan(Info)), sum(Info_pRandSpikePerm<0.01)), 't');
elseif Trill==1
    suplabel(sprintf('%s Coherence %s Trills only session %s N=%d cells %d with spike shuffling test p<0.01', MS, SessionChoice, FeatureName2, sum(~isnan(Info)), sum(Info_pRandSpikePerm<0.01)), 't');
end

fprintf(1,'Cell with highest Info Value: %s', CellsPath(Info == max(Info)).name)
% Amplitude: 59834_20190614_SSS_1-100.mat (cc=91/488) (cc=212/608)
% SpectralMean 59834_20190610_SSS_1-130.mat (cc=43/488)
% Saliency: 59834_20190708_SSM_1-228.mat (cc=454/488)

% Amp Free session: 59834_20190604_SSM_3-229.mat (cc=7/485) (cc = 130/608)
% Amp Operant session: 59834_20190614_SSS_1-100.mat

% Plot the values of the secondary peaks found for some cells
figure(7)
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
    Free = load(fullfile(HDPath,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, 'Free')));
    Operant = load(fullfile(HDPath,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, 'Operant')));
    
elseif Trill==1
    Free = load(fullfile(HDPath,sprintf('%sCoherence_%s_%s_Trill.mat', MS, FeatureName{fn}, 'Free')));
    Operant = load(fullfile(HDPath,sprintf('%sCoherence_%s_%s_Trill.mat', MS, FeatureName{fn}, 'Operant')));
    
end
Colors = [(Operant.Info_pRandSpikePerm<0.01) zeros(length(Operant.Info),1) (Free.Info_pRandSpikePerm<0.01)];
SigCells_logical = ((Operant.Info_pRandSpikePerm<0.01) + (Free.Info_pRandSpikePerm<0.01))>0;
Row59834 = ((Operant.BatID==59834) + (Free.BatID==59834)) >0;
Row11689 = ((Operant.BatID==11689) + (Free.BatID==11689)) >0;

if strcmp(FeatureName{fn}, 'amp')
    FeatureName2 = 'Amplitude';
elseif strcmp(FeatureName{fn}, 'sal')
    FeatureName2 = 'Pitch Saliency';
elseif strcmp(FeatureName{fn}, 'SpectralMean')
    FeatureName2 = 'Spectral Mean';
end
subplot(1,3,1)
scatter(Operant.Info(logical(Row59834)),Free.Info(logical(Row59834)),40,Colors(logical(Row59834) ,:), 'filled')
hold on
scatter(Operant.Info(logical(Row11689)),Free.Info(logical(Row11689)),40,Colors(logical(Row11689) ,:), 'filled','d')
%     scatter(Operant.Info,Free.Info,40,Colors, 'filled')
hold on
plot([0 9], [0 9], 'r:', 'LineWidth',2)
hold off
xlabel('Operant Session')
ylabel('Free Session')
ValidCells4Plot = logical((~isnan(Operant.Info .* Free.Info)));
if ~Trill
    title(sprintf('All Cells (bits/s; N=%d)', sum(ValidCells4Plot)))
    text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
elseif Trill==1
    title(sprintf('All Cells (bits/s; N=%d)', sum(ValidCells4Plot)))
    text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
end

subplot(1,3,2)
scatter(Operant.Info(logical(SigCells_logical.*Row59834)),Free.Info(logical(SigCells_logical.*Row59834)),40,Colors(logical(SigCells_logical.*Row59834) ,:), 'filled')
hold on
scatter(Operant.Info(logical(SigCells_logical.*Row11689)),Free.Info(logical(SigCells_logical.*Row11689)),40,Colors(logical(SigCells_logical.*Row11689) ,:), 'filled','d')
%     scatter(Operant.Info,Free.Info,40,Colors, 'filled')
hold on
plot([0 9], [0 9], 'r:', 'LineWidth',2)
hold off
xlabel('Operant Session')
ylabel('Free Session')
ValidCells4Plot = logical((~isnan(Operant.Info .* Free.Info)).*SigCells_logical);
if ~Trill
    title(sprintf('All Coherence significative Cells (bits/s; N=%d)', sum(ValidCells4Plot)))
    text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
elseif Trill==1
    title(sprintf('All Coherence significative Cells (bits/s; N=%d)', sum(ValidCells4Plot)))
    text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
end


subplot(1,3,3)
scatter(Operant.Info(OpvsFr_Ks2pval_boot>0.01),Free.Info(OpvsFr_Ks2pval_boot>0.01),40,Colors(OpvsFr_Ks2pval_boot>0.01,:), 'filled')
hold on
plot([0 9], [0 9], 'r:', 'LineWidth',2)
hold off
xlabel('Operant Session')
ylabel('Free Session')
ValidCells4Plot = logical((~isnan(Operant.Info .* Free.Info)).*(OpvsFr_Ks2pval_boot>0.01));
if ~Trill
    title(sprintf('Stable Cells with KS test-pvalue>0.01 (bits/s; N=%d)', sum(ValidCells4Plot)))
    text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
    suplabel(sprintf('Information on %s coherence with sound %s of all calls', MS, FeatureName2))
elseif Trill==1
    title(sprintf('Stable Cells with  KStest-pvalue>0.01 (bits/s; N=%d)',sum(ValidCells4Plot)))
    text(0.2, 4.5,  sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))
    suplabel(sprintf('Information on %s coherence with sound %s of Trills only', MS, FeatureName2))
end



figure(8)
clf
[h,pstable]=ttest(Operant.Info(OpvsFr_Ks2pval_boot>0.01), Free.Info(OpvsFr_Ks2pval_boot>0.01))
% [h,pall]=ttest(Operant.Info(SigCells_logical), Free.Info(SigCells_logical))
[h,pall]=ttest(Operant.Info, Free.Info)
% [h,pstable]=kstest(Operant.Info(KDE_OpvsFr_p>0.05)-Free.Info(KDE_OpvsFr_p>0.05))
% [h,pall]=kstest(Operant.Info-Free.Info)
subplot(1,2,1)
histogram(Operant.Info - Free.Info,20, 'FaceColor','k')
% histogram(Operant.Info(SigCells_logical) - Free.Info(SigCells_logical),20, 'FaceColor','k')
xlabel('Information on Motor Coherence difference Operant - Free')
%     ylabel('All cells')
% ylabel('Significant Cells for motor MI coherence during operant OR Free session')
ylabel('# Cells')
xlim([-8 8])
hold on
VL = vline(0, 'r--');
VL.LineWidth=2;
N=sum(~isnan(Operant.Info .* Free.Info));
% N=sum(~isnan(Operant.Info(SigCells_logical) .* Free.Info(SigCells_logical)));
title(sprintf('Information change mean=%.2f+/-%.2f t-test p=%.2f N=%d',nanmean(Operant.Info - Free.Info),nanstd(Operant.Info - Free.Info)./(N)^0.5, pall, N))
% title(sprintf('Information change mean=%.2f+/-%.2f t-test p=%.2f N=%d',nanmean(Operant.Info(SigCells_logical) - Free.Info(SigCells_logical)),nanstd(Operant.Info(SigCells_logical) - Free.Info(SigCells_logical))./(N)^0.5, pall,N ))
subplot(1,2,2)
histogram(Operant.Info(OpvsFr_Ks2pval_boot>0.01) - Free.Info(OpvsFr_Ks2pval_boot>0.01),20, 'FaceColor','k')
xlim([-5 5])
hold on
VL = vline(0, 'r--');
VL.LineWidth=2;
xlabel('Information on Motor Coherence difference Operant - Free')
ylabel('# Stable cells')
NStable =sum((~isnan(Operant.Info .* Free.Info)).*(OpvsFr_Ks2pval_boot>0.01));
title(sprintf('Information change mean=%.2f+/-%.2f t-test p=%.2e N=%d',nanmean(Operant.Info(OpvsFr_Ks2pval_boot>0.01) - Free.Info(OpvsFr_Ks2pval_boot>0.01)),nanstd(Operant.Info(OpvsFr_Ks2pval_boot>0.01) - Free.Info(OpvsFr_Ks2pval_boot>0.01))./(NStable)^0.5, pstable, NStable))
if Trill
    suplabel('Information on coherence of sound Amp with spike rate for Trills only')
else
    suplabel('Information on coherence of sound Amp with spike rate for all calls only')
end



if ~Trill
    MS= 'Auditory';
    Free = load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, 'Free')));
    Operant = load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', MS, FeatureName{fn}, 'Operant')));
    Colors = [(Operant.Info_pTime<0.01) zeros(length(Operant.Info),1) (Free.Info_pTime<0.01)];
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
    text(0.2, 4.5, sprintf('Sig: Red=Operant (N=%d) Blue=Free (N=%d) Magenta=Both (N=%d)', sum(Operant.Info_pTime(ValidCells4Plot)<0.01), sum(Free.Info_pTime(ValidCells4Plot)<0.01),sum(sum(Colors(ValidCells4Plot,:),2)==2)))
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
Auditory = load(fullfile(HDPath,sprintf('AuditoryCoherence_%s_%s.mat', FeatureName{fn}, SessionChoice)));
Motor = load(fullfile(HDPath,sprintf('MotorCoherence_%s_%s.mat', FeatureName{fn}, SessionChoice)));
Row59834 = logical((Auditory.BatID==59834) .* (Motor.BatID==59834));
Row11689 = logical((Auditory.BatID==11689) .* (Motor.BatID==11689));
Colors = [(Motor.Info_pRandSpikePerm<0.01) zeros(length(Motor.Info),1) (Auditory.Info_pRandSpikePerm<0.01)];
if strcmp(FeatureName{fn}, 'amp')
    FeatureName2 = 'Amplitude';
elseif strcmp(FeatureName{fn}, 'sal')
    FeatureName2 = 'Pitch Saliency';
elseif strcmp(FeatureName{fn}, 'SpectralMean')
    FeatureName2 = 'Spectral Mean';
end
subplot(3,1,1)
scatter(Motor.Info(Row59834),Auditory.Info(Row59834),40,Colors(Row59834,:), 'filled')
hold on
scatter(Motor.Info(Row11689),Auditory.Info(Row11689),40,Colors(Row11689, :), 'filled','d')
hold on
plot([0 9], [0 9], 'r:', 'LineWidth',2)
hold off
xlabel('Motor')
ylabel('Auditory')
ValidCells4Plot = ~isnan(Motor.Info .* Auditory.Info);
title(sprintf('Information on Coherence with sound %s during %s Session (bits; N=%d)', FeatureName2, SessionChoice, sum(ValidCells4Plot)))
text(0.2, 7,  sprintf('Sig: Red=Motor (N=%d) Blue=Auditory (N=%d) Magenta=Both (N=%d)', sum(Colors(ValidCells4Plot,1)),sum(Colors(ValidCells4Plot,3)), sum(sum(Colors(ValidCells4Plot,:),2)==2)))

subplot(3,1,2)
[h,pMA]=ttest(Motor.Info, Auditory.Info)
histogram(Motor.Info - Auditory.Info,20, 'FaceColor','k')
xlabel('Information on Coherence difference Motor - Auditory')
ylabel('All cells')
%     xlim([-1.5 1.5])
hold on
VL = vline(0, 'r--');
VL.LineWidth=2;
title(sprintf('Information change in Free session t-test p=%.2e N=%d', pMA, sum(~isnan(Motor.Info .* Auditory.Info))))
%     set(gca, 'Color', 'none')
%     set(gcf, 'Color', 'none')

ErrorMotor = Motor.Info_up - Motor.Info_low;
ErrorAuditory = Auditory.Info_up - Auditory.Info_low;
DPrime = (Motor.Info - Auditory.Info)./((ErrorMotor.^2 + ErrorAuditory.^2).^0.5);
subplot(3,1,3)
histogram(DPrime ,20, 'FaceColor','k')
xlabel('D-prime Information on Coherence difference Motor - Auditory')
ylabel('All cells')
%     xlim([-1.5 1.5])
hold on
VL = vline(0, 'r--');
VL.LineWidth=2;
title(sprintf('D-prime Information change in Free session mean = %.2f +/- %.2f N=%d', nanmean(DPrime),nanstd(DPrime), sum(~isnan(Motor.Info .* Auditory.Info))))

% end
%% effect of sampling on values of Coherence MI
Operant = load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', 'Motor', 'amp', 'Operant')));
figure(); clf;
subplot(1,2,1)
scatter(Operant.Info, Operant.NStims,20,[(Operant.Info_pRandSpikePerm<0.01) zeros(length(Operant.Info),1) zeros(length(Operant.Info),1)], 'filled');
ylabel('# vocalizations')
xlabel('Motor Cohernce MI Operant')
subplot(1,2,2)
scatter(Operant.Info, Operant.LengthX,20,[(Operant.Info_pRandSpikePerm<0.01) zeros(length(Operant.Info),1) zeros(length(Operant.Info),1)], 'filled');
ylabel('Sample size')
xlabel('Motor Cohernce MI Operant')
suplabel(sprintf('Operant Session N= %d/%d significant units', sum(Operant.Info_pRandSpikePerm<0.01), sum(~isnan(Operant.Info_pRandSpikePerm))))


Free = load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', 'Motor', 'amp', 'Free')));
figure(); clf;
subplot(1,2,1)
scatter(Free.Info, Free.NStims,20,[(Free.Info_pRandSpikePerm<0.01) zeros(length(Free.Info),1) zeros(length(Free.Info),1)], 'filled');
ylabel('# vocalizations')
xlabel('Motor Cohernce MI Free')
subplot(1,2,2)
scatter(Free.Info, Free.LengthX,20,[(Free.Info_pRandSpikePerm<0.01) zeros(length(Free.Info),1) zeros(length(Free.Info),1)], 'filled');
ylabel('Sample size')
xlabel('Motor Cohernce MI Free')
suplabel(sprintf('Free Session N= %d/%d significant units', sum(Free.Info_pRandSpikePerm<0.01), sum(~isnan(Free.Info_pRandSpikePerm))))

Free = load(fullfile(Path,sprintf('%sCoherence_%s_%s.mat', 'Auditory', 'amp', 'Free')));
figure(); clf;
subplot(1,2,1)
scatter(Free.Info, Free.NStims,20,[zeros(length(Free.Info),1) (Free.Info_pRandSpikePerm<0.01) zeros(length(Free.Info),1)], 'filled');
ylabel('# vocalizations')
xlabel('Auditory Cohernce MI Free')
subplot(1,2,2)
scatter(Free.Info, Free.LengthX,20,[zeros(length(Free.Info),1) (Free.Info_pRandSpikePerm<0.01) zeros(length(Free.Info),1)], 'filled');
ylabel('Sample size')
xlabel('Auditory Cohernce MI Free')
suplabel(sprintf('Free Session N= %d/%d significant units', sum(Free.Info_pRandSpikePerm<0.01), sum(~isnan(Free.Info_pRandSpikePerm))))

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
%% Explore the neural noise: Poisson or Gamma?!
Session = 'Operant';
Self = 1;
Delay = 200;
FilterSize = 250; % (in ms)
MC = load(fullfile(HDPath,sprintf('MotorCoherence_%s_%s.mat', 'amp', Session)));
GoodInfo=find(~isnan(MC.Info));
% GoodInfo=find(MC.Info_pRandSpikePerm<0.05);
NCells = length(GoodInfo);
CoeffVar_Original = nan(NCells,1);
FanoFactor_Original = cell(NCells,1);
GammaA_Original = nan(NCells,1);
SSEGam_Original = nan(NCells,1);
SSEExp_Original = nan(NCells,1);
CoeffVar_Rescaled = nan(NCells,1);
FanoFactor_Rescaled = cell(NCells,1);
GammaA_Rescaled = nan(NCells,1);
SSEGam_Rescaled = nan(NCells,1);
SSEExp_Rescaled = nan(NCells,1);
SpikeCountPerStim = nan(NCells,2);
GammaA_GLM = nan(NCells,1);
GammaAmpDispersion = nan(NCells,1);
%% The loop
for nc=9:NCells %(demo cell= 560)
    cc=GoodInfo(nc);
%     cc=574;
    fprintf(1,'Cell %d/%d\n',nc,NCells)
    % load data
    Cell = load(fullfile(MC.CellsPath(cc).folder,MC.CellsPath(cc).name));
    if MC.Info(cc)~=Cell.MotorCoherenceOperant.Info
        keyboard
    end
    
    %% Number of vocalizations in the dataset
    if ~isfield(Cell, 'What')
        fprintf(1,'*** . Problem with Cell %d, no what field!! ****\n', cc)
        keyboard
        %         continue
    else
        
        if Self % Select vocalizations from the requested session, emitted by self and that have the right delay before and after
            IndVoc = find(contains(Cell.What, 'Voc') .* contains(Cell.ExpType, Session(1)) .* contains(Cell.Who, 'self') .* (Cell.DelayBefore>=Delay) .* (Cell.DelayAfter>=Delay));
        else % Select vocalizations from the requested session, emitted by...
            % Others with no overlap with other vocalizations and that have
            % the right delay before and after, % ensure as well that there is
            % no noise on the microphone track (manual annotation) if we are
            % in free session, no check for operant
            IndVoc = find(contains(Cell.What, 'Voc') .* contains(Cell.ExpType, Session(1)) .*(~contains(Cell.Who, 'self')) .* (Cell.VocOverlap==0));
            if strcmp(Session(1), 'F')
                IndVoc = intersect(IndVoc, find(Cell.AudioQuality==1));
            end
        end
        NStims = length(IndVoc);
        StimDura = nan(NStims,1);
        for ss=1:NStims
            StimDura(ss) = round(length(Cell.BioSound{IndVoc(ss),2}.sound) ./(Cell.BioSound{IndVoc(ss),2}.samprate)*10^3);
        end
        if any(abs(StimDura - round(Cell.Duration(IndVoc)))>1)
            fprintf(1,'*** . Problem with Cell %d, duration inconcistency!! ****\n', cc)
            keyboard
            %             continue
        end
    end
    
    % Get the optimal Time resolution value given by coherence for
    % Amplitude and spike rate signal
    TR = round(1000./(2*MC.CoherenceWeightedFreq(cc)));
    
    % Compute neural vectors
    % Neural Data loop
    % neural response takes into account all spikes starting
    % at 200 - FilterSize/2 ms before stim onset and stop at 200 - FilterSize/2 ms after
    % stim offset to make sure we're not looking at spike for which we
    % cannot establish the corresponding time varying acoustic feature from
    % -FilterSize/2 ms to FilterSize/2 ms.
    TLim = [FilterSize/2-200 200-FilterSize/2];
    [YPerStim, YPerStimt, YPatterns, SATPerStim, ISIPerStim] = get_y_4GammaFit(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),TLim,TR);
    
    % These are the spike arrival times for each vocalization
    SAT = cellfun(@find, YPatterns, 'UniformOutput', false);
    % Let's keep track of the average number of spikes per stim
    SpikeCountPerStim(nc,1) = nanmean(cellfun(@length, SAT));
    SpikeCountPerStim(nc,2) = nanstd(cellfun(@length, SAT));
    
    % Correct the time of spikes by the time varying rate to see what is
    % the underlying homogenous Gamma/Poisson process
    [ISIPerStim_Rescaled, ISIvecold, YPatterns_Rescaled]=time_rescaled_ISI(YPatterns, YPerStim);
    
    % ISIPerStim are the original inter-spike interval for each vocalization
    if isempty([ISIPerStim{:}]) || length([ISIPerStim{:}])<10
        warning('too few spikes! cannot calculate ISI!/n')
        continue
    end
    % This is the coefficient of variation of ISI
    CoeffVar_Original(nc) = std([ISIPerStim{:}])/mean([ISIPerStim{:}]);
    CoeffVar_Rescaled(nc) = std([ISIPerStim_Rescaled{:}])/mean([ISIPerStim_Rescaled{:}]);
    % This is the Fano Factor of the spike count per bin TR (=1 for a poisson process)
    SpikeCountPerBin_Rescaled = cell(size(YPatterns_Rescaled));
    for ss=1:length(YPatterns_Rescaled)
        VocDuration = length(YPatterns_Rescaled{ss});
        TimeBins = 1:TR:VocDuration;
        SpikeCountPerBin_Rescaled{ss} = nan(length(TimeBins)-1,1);
        for tt=1:(length(TimeBins)-1)
            SpikeCountPerBin_Rescaled{ss}(tt) = sum(YPatterns_Rescaled{ss}(TimeBins(tt):(TimeBins(tt+1)-1)));
        end
    end
    SpikeCountPerBin = cell(size(YPatterns));
    for ss=1:length(YPatterns)
        VocDuration = length(YPatterns{ss});
        TimeBins = 1:TR:VocDuration;
        SpikeCountPerBin{ss} = nan(length(TimeBins)-1,1);
        for tt=1:(length(TimeBins)-1)
            SpikeCountPerBin{ss}(tt) = sum(YPatterns{ss}(TimeBins(tt):(TimeBins(tt+1)-1)));
        end
    end
    FanoFactor_Original{nc} = ((cellfun(@std, SpikeCountPerBin)).^2) ./ cellfun(@mean, SpikeCountPerBin);
    FanoFactor_Rescaled{nc} = ((cellfun(@std, SpikeCountPerBin_Rescaled)).^2) ./ cellfun(@mean, SpikeCountPerBin_Rescaled);
    % plot the histogram of ISI
    figure(20);subplot(2,1,1);cla
    H=histogram([ISIPerStim{:}], 'BinWidth',1); xlabel('ISI (ms)');
    % fit the distribution of ISI with a gamma distribution
    Gam = fitdist([ISIPerStim{:}]', 'Gamma');
    GammaA_Original(nc) = Gam.a;
    % fit the distribution of ISI with an exponential (Poisson process)
    Poi = fitdist([ISIPerStim{:}]', 'Exponential');
    % add the2 fit to the figure
    XISI = H.BinEdges(1):(H.BinEdges(end)-1);
    Xexp = pdf('Exponential',XISI,Poi.mu);
    Xgam = gampdf(XISI, Gam.a, Gam.b);
    % Calculate SSE
    PdfH = H.BinCounts./sum(H.BinCounts);
    SSEGam_Original(nc) = sum((PdfH-Xgam).^2);
    SSEExp_Original(nc) = sum((PdfH-Xexp).^2);
    % plot everything
    figure(20);subplot(2,1,1); cla
    yyaxis right; plot(XISI, Xexp, '-', 'LineWidth',2, 'Color', [0.929 0.694 0.125]);
    hold on; yyaxis right; plot(XISI, Xgam, '-', 'LineWidth',2);
    hold on; yyaxis left; H=histogram([ISIPerStim{:}], 'BinWidth',1); xlabel('Original ISI (ms)');hold off
    legend({'ISI','Exp fit', 'Gamma fit'})
    title(sprintf('ISI fit Gamma parameter = %.2f SSEGam = %.2e SSEExp = %.2e FanoFactor = %.2f +/- %.2f', Gam.a, SSEGam_Original(nc), SSEExp_Original(nc), nanmean(FanoFactor_Original{nc}), nanstd(FanoFactor_Original{nc})/(length(FanoFactor_Original{nc}))^0.5))
    suplabel('Original data', 't')
    
    % plot the histogram of ISI Rescaled
    figure(20);subplot(2,1,2);cla
    H=histogram([ISIPerStim_Rescaled{:}], 'BinWidth',1); xlabel('ISI Rescaled (ms)');
    % fit the distribution of ISI with a gamma distribution
    Gam = fitdist([ISIPerStim_Rescaled{:}]', 'Gamma');
    GammaA_Rescaled(nc) = Gam.a;
    % fit the distribution of ISI with an exponential (Poisson process)
    Poi = fitdist([ISIPerStim_Rescaled{:}]', 'Exponential');
    % add the2 fit to the figure
    XISI = H.BinEdges(1):(H.BinEdges(end)-1);
    Xexp = pdf('Exponential',XISI,Poi.mu);
    Xgam = gampdf(XISI, Gam.a, Gam.b);
    % Calculate SSE
    PdfH = H.BinCounts./sum(H.BinCounts);
    SSEGam_Rescaled(nc) = sum((PdfH-Xgam).^2);
    SSEExp_Rescaled(nc) = sum((PdfH-Xexp).^2);
    figure(20);subplot(2,1,2);cla
    yyaxis right; plot(XISI, Xexp, '-', 'LineWidth',2, 'Color', [0.929 0.694 0.125]);
    hold on; yyaxis right; plot(XISI, Xgam, '-', 'LineWidth',2);
    hold on; yyaxis left; H=histogram([ISIPerStim_Rescaled{:}], 'BinWidth',1); xlabel('ISI Rescaled (ms)');hold off
    legend({'ISI Rescaled','Exp fit', 'Gamma fit'})
    title(sprintf('ISI Rescaled fit Gamma parameter = %.2f SSEGam = %.2e SSEExp = %.2e FanoFactor = %.2f +/- %.2f', Gam.a, SSEGam_Rescaled(nc), SSEExp_Rescaled(nc), nanmean(FanoFactor_Rescaled{nc}), nanstd(FanoFactor_Rescaled{nc})/(length(FanoFactor_Rescaled{nc}))^0.5))
    suplabel(sprintf('Cell ID: %s', MC.CellsPath(cc).name), 't')
    
    % Estimate the shape parameter of gamma distribution based on a linear
    % model with Amplitude 
     % Calculate acoustic features input to the models
        % acoustic data is a Matrix of the value of the acoustic feature sampled
        % at 1000Hz that for each Spike (at t) starts
        % t-FilterSize/2 and stops at
        % t+FilterSize/2 with t references to zero at stim onset
        DefaultVal = 0;%zero should be the default value for the amplitude, we know here that there is no sound
        [XAmpPerStim]  = get_x_4gammafit(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), SATPerStim, FilterSize,DefaultVal,'amp');
        XAmp = [XAmpPerStim{:}]';
            
        [~, XAmp_Mu, XAmp_Sigma] = zscore(reshape(XAmp,numel(XAmp),1));
        XAmp_ZS = (XAmp - XAmp_Mu) ./ XAmp_Sigma;
        if any(isnan(XAmp_ZS))
            keyboard
        end
        [AmpPC,AmpScore,AmpLatent,~,AmpPCAExplained,AmpPCAMu] = pca(XAmp_ZS);
        NPC = find(cumsum(AmpPCAExplained)>95,1);
        % Run the gamma model
        % Now let's fit the data with a GLM gamma
        [B_GLM,Dev, Stats] = glmfit(AmpScore(:,1:NPC),[ISIPerStim{:}],'gamma', 'link', 'log');
        GamModel = fitglm(AmpScore(:,1:NPC),[ISIPerStim{:}],'linear','Distribution','gamma', 'DispersionFlag', true,'Link', 'log');
        
        % estimated shape
        fprintf(1,'The model estimates a shape of %.2f\n', 1/GamModel.Dispersion)
        fprintf(1, 'The model estimates the intercept to be %.1f\n', GamModel.Coefficients.Estimate(1) - log(1/GamModel.Dispersion))
        fprintf(1, 'The model estimates the Beta coefficients to be:\n')
        GammaA_GLM(nc) = 1/GamModel.Dispersion;
        GammaAmpDispersion(nc) = GamModel.Dispersion;
        GamModel.Coefficients.Estimate
        figure(20)
        suplabel(sprintf('GLM Gamma shape = %.2f', GammaA_GLM(nc)))
        drawnow
%     pause()
end
save(fullfile(HDPath, ('OperantVocPoissonGammaNeuralNoise.mat')), 'CoeffVar_Original', 'FanoFactor_Original', 'GammaA_Original', 'SSEGam_Original', 'SSEExp_Original', 'CoeffVar_Rescaled', 'FanoFactor_Rescaled', 'GammaA_Rescaled', 'SSEGam_Rescaled', 'SSEExp_Rescaled', 'SpikeCountPerStim','GammaA_GLM', 'GammaAmpDispersion')

%% Plot figures of Poisson vs Gamma cell behavior
Session = 'Operant';
MC = load(fullfile(HDPath,sprintf('MotorCoherence_%s_%s.mat', 'amp', Session)));
NPoissonCells_Original = sum((cellfun(@nanmean, FanoFactor_Original)<1.02) .* (cellfun(@nanmean, FanoFactor_Original)>0.98) .* (CoeffVar_Original>0.95) .* (CoeffVar_Original<1.05));
NTotCells = sum(~isnan(CoeffVar_Original));
FanoFactorAv_Original = cellfun(@nanmean, FanoFactor_Original);
figure(21);clf; subplot(2,3,1); cubehelix_niceplot(cellfun(@nanmean, FanoFactor_Original), GammaA_Original, SpikeCountPerStim(:,1),3, ' '); xlabel('FanoFactor (var/mean)'); ylabel('Gamma fit shape');
subplot(2,3,2); cubehelix_niceplot(CoeffVar_Original, GammaA_Original, SpikeCountPerStim(:,1),3, ' ' ); xlabel('Coefficient of variation std/mean'); ylabel('Gamma fit shape');
subplot(2,3,4); cubehelix_niceplot(FanoFactorAv_Original, CoeffVar_Original, SpikeCountPerStim(:,1),3, sprintf('%d cells are Poisson like (%f.2 %%)', NPoissonCells_Original, NPoissonCells_Original*100/NTotCells)); xlabel('FanoFactor (var/mean)'); ylabel('Coefficient of variation std/mean'); hold on; vline([0.98 1.02], {'k', 'k'});hold on;  hline([0.95 1.05], {'k', 'k'})
subplot(2,3,5); cubehelix_niceplot(FanoFactorAv_Original(~isnan(GammaA_Original)), CoeffVar_Original(~isnan(GammaA_Original)), GammaA_Original(~isnan(GammaA_Original)),3, sprintf('%d cells are Poisson like (%f.2 %%)', NPoissonCells_Original, NPoissonCells_Original*100/NTotCells)); xlabel('FanoFactor (var/mean)'); ylabel('Coefficient of variation std/mean'); hold on; vline([0.98 1.02], {'k', 'k'});hold on;  hline([0.95 1.05], {'k', 'k'})
subplot(2,3,6); cubehelix_niceplot(MC.Info(GoodInfo), CoeffVar_Original, SpikeCountPerStim(:,1),3); xlabel(' Info on Motor coherence with amplitude');ylabel('Coefficient of variation std/mean');
suplabel('Original data','t')
FanoFactorAv_Rescaled = cellfun(@nanmean, FanoFactor_Rescaled);
NPoissonCells_Rescaled = sum((cellfun(@nanmean, FanoFactor_Rescaled)<1.02) .* (cellfun(@nanmean, FanoFactor_Rescaled)>0.98) .* (CoeffVar_Rescaled>0.95) .* (CoeffVar_Rescaled<1.05));
figure(22);clf; subplot(2,3,1); cubehelix_niceplot(cellfun(@nanmean, FanoFactor_Rescaled), GammaA_Rescaled, SpikeCountPerStim(:,1),3, ' '); xlabel('FanoFactor (var/mean)'); ylabel('Gamma fit shape');
subplot(2,3,2); cubehelix_niceplot(CoeffVar_Rescaled, GammaA_Rescaled, SpikeCountPerStim(:,1),3, ' ' ); xlabel('Coefficient of variation std/mean'); ylabel('Gamma fit shape');
subplot(2,3,4); cubehelix_niceplot(FanoFactorAv_Rescaled, CoeffVar_Rescaled, SpikeCountPerStim(:,1),3, sprintf('%d cells are Poisson like (%f.2 %%)', NPoissonCells_Rescaled, NPoissonCells_Rescaled*100/NTotCells)); xlabel('FanoFactor (var/mean)'); ylabel('Coefficient of variation std/mean'); hold on; vline([0.98 1.02], {'k', 'k'});hold on;  hline([0.95 1.05], {'k', 'k'})
subplot(2,3,5); cubehelix_niceplot(FanoFactorAv_Rescaled(~isnan(GammaA_Rescaled)), CoeffVar_Rescaled(~isnan(GammaA_Rescaled)), GammaA_Rescaled(~isnan(GammaA_Rescaled)),3, sprintf('%d cells are Poisson like (%f.2 %%)', NPoissonCells_Rescaled, NPoissonCells_Rescaled*100/NTotCells)); xlabel('FanoFactor (var/mean)'); ylabel('Coefficient of variation std/mean'); hold on; vline([0.98 1.02], {'k', 'k'});hold on;  hline([0.95 1.05], {'k', 'k'})
subplot(2,3,6); cubehelix_niceplot(MC.Info(GoodInfo), CoeffVar_Rescaled, SpikeCountPerStim(:,1),3); xlabel(' Info on Motor coherence with amplitude');ylabel('Coefficient of variation std/mean');
suplabel('Rescaled data', 't')

figure(); subplot(1,3,1);scatter(CoeffVar_Original, CoeffVar_Rescaled, 'filled'); hold on ; line([0 2], [0 2], 'Color','k');xlabel('Coefficient of variation original'); ylabel('Coeffificient of variation rescaled')
subplot(1,3,2);scatter(cellfun(@nanmean,FanoFactor_Original), cellfun(@nanmean,FanoFactor_Rescaled), 'filled'); hold on ; line([0.93 1.05], [0.93 1.05], 'Color','k');xlabel('FanoFactor original'); ylabel('FanoFactor rescaled')
subplot(1,3,3);scatter(GammaA_Original, GammaA_Rescaled, 'filled'); hold on ; line([0 45], [0 45], 'Color','k');xlabel('GammaA original'); ylabel('GammaA rescaled')


%% Gamma model for all cells
% now we have a better estimate of the shape for the Gamma noise
        % (shape obtained with the GLM fit)
        % Let's use that shape to get a better estimate for the time
        % varying spike rate (let's say that we fix the CV of the estimate
        % of the ISI across time (CVe = SE/mean = 0.1 or 0.5?), then given
        % the shape of the gamma distribution of the noise, how many spikes
        % do we need in nearest neighbour spike rate estimate to get the
        % bet estimate of the rate. Because the mean of the distribution of
        % ISI is the same as the mean of your estimates of the mean, and
        % mean = shape * scale at any time point and var = shape * scale^2,
        % then SE = (var/n)^0.5 = scale * (shape/n)^0.5 = mean / (shape *
        % n)^0.5 -> n = (mean/SE)^2 * 1/shape = 1/(shape * CVe^2)
        CVe = 0.2; % This is the coeffcient of variation of the estimate of the mean interspike interval, we constrain the system by deciding to fix it at 0.2 (that's the noise on the calculatoin of the mean) 
        NumNearNeighSpike = round(GamModel.Dispersion/(CVe.^2));
        [YPerStim, YPerStimt, YPatterns, SATPerStim, ISIPerStim] = get_y_4GammaFit(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),TLim,TR, NumNearNeighSpike);
        
        % Now that we have a better estimate of the time varying rate, we
        % need to generate spike patterns with gamma noise that follow that
        % rate for each vocalization


%% Explore Cell by cell the profile of coherence and scatter plots of cell tuning for Good Cells
%% Then for all cells
InputFig=0;
Session = 'Operant';
Self = 1;
FilterSize = 250; % (in ms)
Delay = 200;
Response_samprate = 1000; % sampling rate of neural data in Hz
MC = load(fullfile(HDPath,sprintf('MotorCoherence_%s_%s.mat', 'amp', Session)));
GoodInfo=find(~isnan(MC.Info));
% GoodInfo=find(MC.Info_pRandSpikePerm<0.05);
NCells = length(GoodInfo);
AmpPredictor = nan(NCells,3); % first element is the F statistic with full model, second is pvalue, third is the partial adjusted R2
SalPredictor = nan(NCells,3);% first element is the F statistic with full model, second is pvalue, third is the partial adjusted R2
SpecMeanPredictor = nan(NCells,3); % first element is the F statistic with full model, second is pvalue, third is the partial adjusted R2
CallTypePredictor = nan(NCells,3); % first element is the F statistic with full model, second is pvalue, third is the partial adjusted R2
FullModelR2 = nan(NCells,4); % first element is the F statistic with Null model, second is pvalue, third is the partial adjusted R2
ModelAmpR2 = nan(NCells,3); % first element is the F statistic with Null model, second is pvalue, third is the adjusted R2
ModelAmpLM_Coeff = cell(NCells,2); % first column are time point X-Y of each coefficient (filter) that are reported in column 2 for each cell for the model with Amp
ModelAmpGLM_Coeff = cell(NCells,2);% first column are time point X-Y of each coefficient (filter) that are reported in column 2 for each cell for the model with Amp
ModelAmpR2_cv = cell(NCells,1); % for each cell, vector of R2 distribution obtain in leave one vocalization out cross-validation with lm
ModelGLMAmpDev_cv = cell(NCells,1); % for each cell, vector of Deviance distribution obtain in leave one vocalization out cross-validation with GLM
ModelGLMAmpGoodnessOfFit = nan(NCells,2); % Goodness of fit of the GLM Poisson model with Amplitude obtained on average in cross-validation (leave one vocalization out; column 1) with its STE (column 2)
ModelCeilR2_cv = cell(NCells,1); % for each cell, vector of R2 distribution obtain in leave one vocalization out cross-validation with lm with ceiling model
ModelGLMCeilDev_cv = cell(NCells,1); % for each cell, vector of Deviance distribution obtain in leave one vocalization out cross-validation with GLM with ceiling model
ModelGLMFloorDev_cv = cell(NCells,1); % for each cell, vector of Deviance distribution obtain in leave one vocalization out cross-validation with GLM with Floor model = average spike count

NBoots = 100;
R2CeilBoot = nan(NCells,NBoots);
DevianceCeilBoot = nan(NCells,NBoots);
SSECeil = nan(NCells,NBoots);
ModelSalR2 = nan(NCells,4); % first element is the F statistic with Null model, second is pvalue, third is the adjusted R2, Fourth is the coefficient in the model
ModelSpecMeanR2 = nan(NCells,4); % first element is the F statistic with Null model, second is pvalue, third is the adjusted R2, Fourth is the coefficient in the model
ModelCallTypeR2 = nan(NCells,4); % first element is the F statistic with Null model, second is pvalue, third is the adjusted R2, Fourth is the coefficient in the model
BestShift = nan(NCells,1); % Optimal delay used in the model to center the filter
%% The loop
for nc=560:NCells
    cc=GoodInfo(nc);
%     cc=574;
    fprintf(1,'Cell %d/%d\n',nc,NCells)
    
    % Plot the scatter tuning curves
    % load data
    Cell = load(fullfile(MC.CellsPath(cc).folder,MC.CellsPath(cc).name));
    if MC.Info(cc)~=Cell.MotorCoherenceOperant.Info
        keyboard
    end
    
    %% Number of vocalizations in the dataset
    if ~isfield(Cell, 'What')
        fprintf(1,'*** . Problem with Cell %d, no what field!! ****\n', cc)
        keyboard
%         continue
    else
        
        if Self % Select vocalizations from the requested session, emitted by self and that have the right delay before and after
            IndVoc = find(contains(Cell.What, 'Voc') .* contains(Cell.ExpType, Session(1)) .* contains(Cell.Who, 'self') .* (Cell.DelayBefore>=Delay) .* (Cell.DelayAfter>=Delay));
        else % Select vocalizations from the requested session, emitted by...
            % Others with no overlap with other vocalizations and that have
            % the right delay before and after, % ensure as well that there is
            % no noise on the microphone track (manual annotation) if we are
            % in free session, no check for operant
            IndVoc = find(contains(Cell.What, 'Voc') .* contains(Cell.ExpType, Session(1)) .*(~contains(Cell.Who, 'self')) .* (Cell.VocOverlap==0));
            if strcmp(Session(1), 'F')
                IndVoc = intersect(IndVoc, find(Cell.AudioQuality==1));
            end
        end
        NStims = length(IndVoc);
        StimDura = nan(NStims,1);
        for ss=1:NStims
            StimDura(ss) = round(length(Cell.BioSound{IndVoc(ss),2}.sound) ./(Cell.BioSound{IndVoc(ss),2}.samprate)*10^3);
        end
        if any(abs(StimDura - round(Cell.Duration(IndVoc)))>1)
            fprintf(1,'*** . Problem with Cell %d, duration inconcistency!! ****\n', cc)
            keyboard
%             continue
        end
    end
    
    % Get the type of the vocalization
    Trill1_Bark0 = contains(Cell.What(IndVoc), 'Tr');
    
    % Get the optimal Time resolution value
    TR = round(1000./(2*MC.CoherenceWeightedFreq(cc)));
    NBins = floor(FilterSize/(2*TR)); % Number of time bins before and after optimal Delay (CoherencyT_DelayAtzero)
    % Compute neural vectors
    % Neural Data loop
    % neural response is a vector that compile all spike counts starting
    % at 200ms (-Delay(1)) before stim onset and stop at 200ms (Delay(2)) after
    % stim offset
        
        if isnan(MC.CoherencyT_DelayAtzero(cc)) % Let's try to figure out a value for this
            figure(6);clf; plot(Cell.MotorCoherenceOperant.CoherencyT_xTimeDelay, Cell.MotorCoherenceOperant.CoherencyT_filt, 'LineWidth',2);xlabel('Time Delay');ylabel('Coherency')
            [P,Locs] = findpeaks(-Cell.MotorCoherenceOperant.CoherencyT_filt);
            
            if ~isempty(Locs)
                [P,IndM] = max(P);
                Locs = Locs(IndM);
                CoherencyT_DelayAtzero = Cell.MotorCoherenceOperant.CoherencyT_xTimeDelay(Locs);
                CrossZero1 = Cell.MotorCoherenceOperant.CoherencyT_xTimeDelay(find(Cell.MotorCoherenceOperant.CoherencyT_filt(1:Locs)<=mean(Cell.MotorCoherenceOperant.CoherencyT_filt), 1, 'last')+1);
                CrossZero2 = Cell.MotorCoherenceOperant.CoherencyT_xTimeDelay(Locs + find(Cell.MotorCoherenceOperant.CoherencyT_filt(Locs+1:end)<=mean(Cell.MotorCoherenceOperant.CoherencyT_filt), 1, 'first')-1);
                if isempty(CrossZero2)
                    CrossZero2 = Cell.MotorCoherenceOperant.CoherencyT_xTimeDelay(Locs + find(Cell.MotorCoherenceOperant.CoherencyT_filt(Locs+1:end)==min(Cell.MotorCoherenceOperant.CoherencyT_filt(Locs+1:end)), 1, 'first')-1);
                end
                CoherencyT_WidthAtMaxPeak = CrossZero2 - CrossZero1;
                
                warning('No delay was originally found for that cell: hyp: inverse coherency curve?')

            else
                warning('No delay was originally found for that cell: hyp: Very Large frequency?')
                CoherencyT_DelayAtzero=Cell.MotorCoherenceOperant.CoherencyT_xTimeDelay(Cell.MotorCoherenceOperant.CoherencyT_filt == max(Cell.MotorCoherenceOperant.CoherencyT_filt));
                keyboard
%                 continue
            end
            if abs(CoherencyT_DelayAtzero)<=200
                TLim = [CoherencyT_DelayAtzero-200+NBins*TR CoherencyT_DelayAtzero+200-NBins*TR];
                BestShift(cc) = CoherencyT_DelayAtzero;
%               TLim = [-MC.CoherencyT_DelayAtzero(cc) MC.CoherencyT_DelayAtzero(cc)];
            else % this is an abherent delay because there is only 200ms before and after each vocaliztaion when calculating the coherence so center the filter at 0
                TLim = [NBins*TR-200 200-NBins*TR];
                BestShift(cc) = 0;
            end
        else
            if abs(MC.CoherencyT_DelayAtzero(cc))<=200
                TLim = [MC.CoherencyT_DelayAtzero(cc)-200+NBins*TR MC.CoherencyT_DelayAtzero(cc)+200-NBins*TR];
%               TLim = [-MC.CoherencyT_DelayAtzero(cc) MC.CoherencyT_DelayAtzero(cc)];
                BestShift(cc) = MC.CoherencyT_DelayAtzero(cc);
            else % this is an abherent delay because there is only 200ms before and after each vocaliztaion when calculating the coherence so center the filter at 0
                TLim = [NBins*TR-200 200-NBins*TR];
                BestShift(cc) = 0;
            end
        end
        [YPerStim, YPerStimt, YPatterns,YCountPerStim] = get_y_4Coherence(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),TLim,TR);
        Y = [YPerStim{:}]';
        YCount = [YCountPerStim{:}]';
        [Y_ZS, Y_Mu, Y_Sigma] = zscore(Y);
        
        if isempty(Y)
            fprintf(1,'no spike during vocalization! No model!\n')
%             continue
        end
        % Get the vector of call type
        Ind = [0 cumsum(cellfun('length',YPerStim))];
        Trill1_Bark0_local = nan(size(Y));
        for ss = 1:NStims
            Trill1_Bark0_local((Ind(ss)+1):Ind(ss+1)) = Trill1_Bark0(ss) .* ones(Ind(ss+1)-Ind(ss),1);
        end
        
        
        % Calculate acoustic features input to the models
        % acoustic data is a Matrix of the value of the acoustic feature sampled
        % at 1000Hz that for each row (at t) starts
        % t-CoherencyT_DelayAtzero-NBins*TR and stops at
        % t-CoherencyT_DelayAtzero+NBins*TR with t references to zero at stim onset
        DefaultVal = 0;%zero should be the default value for the amplitude, we know here that there is no sound
        XAmp = nan(length(Y), 2*NBins+1);
        Bins = -NBins:NBins;
        for bb = 1:(2*NBins+1)
            if ((Bins(bb))*TR - BestShift(cc))<0
                fprintf(1, 'Calculating the Regressor Amplitude %d/%d that anticipates neural response by %d ms\n', bb, (2*NBins+1),(Bins(bb))*TR- BestShift(cc));
            elseif ((Bins(bb))*TR - BestShift(cc))>0
                fprintf(1, 'Calculating the Regressor Amplitude %d/%d that follows neural response by %d ms\n', bb, (2*NBins+1), (Bins(bb))*TR- BestShift(cc));
            elseif ((Bins(bb))*TR - BestShift(cc))==0
                fprintf(1, 'Calculating the Regressor Amplitude %d/%d that is just alligned at stimulus onset with neural response\n', bb, (2*NBins+1));
            end
            Range = TLim - BestShift(cc) + (Bins(bb))*TR;
            [XAmpPerStim,XAmpPerStimt]  = get_x_4coherence(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Range,TR,DefaultVal,'amp');
            XAmp(:,bb) = [XAmpPerStim{:}]';
            if InputFig
                plotxyfeaturescoherence(Cell.BioSound(IndVoc,2),YPerStim,YPerStimt,XAmpPerStim,XAmpPerStimt,TR,Range,Cell.Duration(IndVoc), 10000,'Amplitude')
            end
        end
        [~, XAmp_Mu, XAmp_Sigma] = zscore(reshape(XAmp,numel(XAmp),1));
        XAmp_ZS = (XAmp - XAmp_Mu) ./ XAmp_Sigma;
        if any(isnan(XAmp_ZS))
            keyboard
        end
        [AmpPC,AmpScore,AmpLatent,~,AmpPCAExplained,AmpPCAMu] = pca(XAmp_ZS);
        NPC = find(cumsum(AmpPCAExplained)>95,1);
        figure(17);clf;plot(cumsum(AmpPCAExplained), 'k-', 'LineWidth',2); ylabel('Explained variance of Amp by PCA'); xlabel('# of Principal Components')
        hold on; vline(NPC, 'g', sprintf('Optimal # PC = %d', NPC));
        
        
        DefaultVal = 0;%zero should be the default value for the saliency, we know here that there is no sound
        [XSalPerStim, XSalPerStimt] = get_x_4coherence(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Range,TR,DefaultVal,'sal', TR/2);
        XSal = [XSalPerStim{:}]';
        [XSal_ZS, XSal_Mu, XSal_Sigma] = zscore(XSal);
        if any(isnan(XSal_ZS))
            keyboard
        end
        
%         DefaultVal = 0;%zero should be the default value for the Fundamental, we know here that there is no sound
%         [XFundPerStim, XFundPerStimt] = get_x_4coherence(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Delay,TR,DefaultVal,'f0', TR/2);
%         XFund = [XFundPerStim{:}]';
%         
        DefaultVal = 'mean';%mean of specmean should be the default value for the specmean, we know here that there is no sound
        [XSpecMeanPerStim, XSpecMeanPerStimt] = get_x_4coherence(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Range,TR,DefaultVal,'SpectralMean', TR/2);
        XSpecMean = [XSpecMeanPerStim{:}]';
        XSpecMean(isnan(XSpecMean)) = nanmean(XSpecMean);
        [XSpecMean_ZS, XSpecMean_Mu, XSpecMean_Sigma] = zscore(XSpecMean);
        
        
        if InputFig
            plotxyfeaturescoherence(Cell.BioSound(IndVoc,2),YPerStim,YPerStimt,XSpecMeanPerStim,XSpecMeanPerStimt,TR,Range,Cell.Duration(IndVoc), 10000,'Spectral Mean')
            plotxyfeaturescoherence(Cell.BioSound(IndVoc,2),YPerStim,YPerStimt,XSalPerStim,XSalPerStimt,TR,Range,Cell.Duration(IndVoc), 10000,'Saliency')
            plotxyfeaturescoherence(Cell.BioSound(IndVoc,2),YPerStim,YPerStimt,XAmpPerStim,XAmpPerStimt,TR,Range,Cell.Duration(IndVoc), 10000,'Amplitude')
%             plotxyfeaturescoherence(Cell.BioSound(IndVoc,2),YPerStim,YPerStimt,XFundPerStim,XFundPerStimt,TR,Range,Cell.Duration(IndVoc), 10000,'Fundamental')
        end
      
        % Run some GLM to establish the effect of acoustic features
        % formula Y~ XAmp + XSal + XSpecMean
        %         TermsMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 0 ];
        TableFull = table(XAmp_ZS, XSal_ZS, XSpecMean_ZS, categorical(Trill1_Bark0_local),Y_ZS,'VariableNames',{'Amplitude' 'Saliency' 'SpectralMean' 'CallType' 'Rate'});
        FullModel = fitlm(TableFull, 'Rate ~ Amplitude + Saliency + SpectralMean + CallType + Amplitude:Saliency + Amplitude:SpectralMean + Amplitude:CallType + Saliency:SpectralMean + Saliency:CallType + SpectralMean:CallType' , 'CategoricalVars', 4);
        FullModelR2(nc,3) = FullModel.Rsquared.Adjusted;
        AnovaF = anova(FullModel, 'summary');
        FullModelR2(nc,1) = AnovaF.F(contains(AnovaF.Properties.RowNames, 'Model'));
        FullModelR2(nc,2) = AnovaF.pValue(contains(AnovaF.Properties.RowNames, 'Model'));
        ModelwoAmp = fitlm(TableFull, 'Rate ~ Saliency + SpectralMean + CallType + Saliency:SpectralMean + Saliency:CallType + SpectralMean:CallType' , 'CategoricalVars', 4);
        
        
        % Run a leave one out vocalization cross validation to get the goodness of fit
        ModelAmpR2_cv{nc} = nan(size(YPerStim));
        ModelGLMAmpDev_cv{nc} = nan(size(YPerStim));
        FirstIndStim = [1 cumsum(cellfun(@length,YPerStim))+1]; FirstIndStim(end) = [];
        LastIndStim = cumsum(cellfun(@length,YPerStim));
        for Stim=1:NStims
            % Linear model approach
            Y_ZS_cv = Y_ZS;
            Y_ZS_cv(FirstIndStim(Stim):LastIndStim(Stim)) = [];
            XAmpScore_cv = AmpScore(:,1:NPC);
            XAmpScore_cv(FirstIndStim(Stim):LastIndStim(Stim),:) = [];
            ModelAmp_cv = fitlm(XAmpScore_cv, Y_ZS_cv );
            Y_ZS_pred = predict(ModelAmp_cv, AmpScore(FirstIndStim(Stim):LastIndStim(Stim),1:NPC));
            SSE0_local = sum((Y_ZS(FirstIndStim(Stim):LastIndStim(Stim)) - mean(Y_ZS_cv)).^2); % The null model here is the average rate accrross all other vocalizations and time points
            ModelAmpR2_cv{nc}(Stim) = 1 - sum((Y_ZS(FirstIndStim(Stim):LastIndStim(Stim)) - Y_ZS_pred).^2)/SSE0_local;
            % Generalized linear model apporach
            YCount_cv = YCount;
            YCount_cv(FirstIndStim(Stim):LastIndStim(Stim)) = [];
            XAmpScore_cv = AmpScore(:,1:NPC);
            XAmpScore_cv(FirstIndStim(Stim):LastIndStim(Stim),:) = [];
            [ModelGLMAmp_cvB,~] = glmfit(XAmpScore_cv, YCount_cv, 'Poisson');
            YCount_pred = glmval(ModelGLMAmp_cvB, AmpScore(FirstIndStim(Stim):LastIndStim(Stim),1:NPC), 'log');
            LL_Saturated = LL_Calculus(YCount(FirstIndStim(Stim):LastIndStim(Stim)),YCount(FirstIndStim(Stim):LastIndStim(Stim)));
            LL_amp = LL_Calculus(YCount(FirstIndStim(Stim):LastIndStim(Stim)), YCount_pred);
            ModelGLMAmpDev_cv{nc}(Stim) = -2*(LL_amp - LL_Saturated);
            % Get the Floor model for the GLM approach
            LL_Floor = LL_Calculus(YCount(FirstIndStim(Stim):LastIndStim(Stim)), round(mean(YCount_cv)));
            ModelGLMFloorDev_cv{nc}(Stim) = -2*(LL_Floor - LL_Saturated);
            % Get the ceiling model value for both approach
            [~, Spike_arrival] = poisson_spiketrain_gen(YPerStim(Stim), YPerStimt(Stim),TLim, TR, Response_samprate, NBoots);
            % Run through the number of bootstrap loops to calculate ceil R2 and ceil deviance
            R2Ceil_local = nan(NBoots,1);
            DevianceCeil_local = nan(NBoots,1);
            DevianceFloor_local = nan(NBoots,1);
            for tt=1:size(Spike_arrival,2)
                [YPerStimbb, YPerStimtbb,~,YCountPerStimbb] = get_y_4Coherence(Spike_arrival(1,tt), Cell.Duration(IndVoc(Stim)),TLim,TR);
                Y_bb = ([YPerStimbb{:}]' - Y_Mu)/Y_Sigma;
                SSECeil_local = sum((Y_ZS(FirstIndStim(Stim):LastIndStim(Stim)) - Y_bb).^2);
                R2Ceil_local(tt) = 1-SSECeil_local/SSE0_local;
                YCountbb = [YCountPerStimbb{:}]';
                LL_Ceilbb = LL_Calculus(YCount(FirstIndStim(Stim):LastIndStim(Stim)), YCountbb);
                DevianceCeil_local(tt) = -2*(LL_Ceilbb - LL_Saturated);
                if InputFig
                    figure(16);
                    clf
                    Win = YPerStimt{Stim} + TR - TLim(1);
                    hist(Spike_arrival{1,tt},Win(end))
                    title(sprintf('NSpike = %d, SSE0 = %.2f SSECeil = %.2f R2 = %.2f DevianceCeil = %.2f',length(Spike_arrival{1,tt}), SSE0_local, SSECeil_local,R2Ceil_local(tt),DevianceCeil_local(tt)))
                    xlim([0 Win(end)]+TLim(1))
                    h=findobj(gca, 'Type', 'patch');
                    h.FaceColor = 'k';
                    h.EdgeColor = 'k';
                    hold on
                    yyaxis right
                    plot(YPerStimt{Stim},YPerStim{Stim}, 'LineWidth',2)
                    ylabel('Original rate mHz')
                    hold on
                    plot(YPerStimtbb{1},YPerStimbb{1}, 'LineWidth',2)
                    hold on
                    Hmean=hline(mean(Y));
                    Hmean.LineWidth = 2;
                    pause
                end
                
            end
            ModelCeilR2_cv{nc}(Stim) = mean(R2Ceil_local);
            ModelGLMCeilDev_cv{nc}(Stim) = mean(DevianceCeil_local);
        end
        figure(6);clf; subplot(2,2,3);histogram(ModelAmpR2_cv{nc}, 'FaceColor', 'k', 'EdgeColor', 'k'); ylabel('# CV folds'); xlabel('Cross validated R2 for LM (leave one vocalization out)'); %Xlim=get(gca, 'XLim'); set(gca, 'XLim', [0 Xlim(2)]);
        hold on; MeanLine = vline(mean(ModelAmpR2_cv{nc})); MeanLine.LineWidth=2; MeanLine.LineStyle = '-';
        hold on;STELine1=vline(mean(ModelAmpR2_cv{nc})-2*std(ModelAmpR2_cv{nc})/NStims^0.5); STELine1.LineWidth=2;
        hold on;STELine2=vline(mean(ModelAmpR2_cv{nc})+2*std(ModelAmpR2_cv{nc})/NStims^0.5);STELine2.LineWidth=2;
        subplot(2,2,4);histogram(ModelGLMAmpDev_cv{nc}, 'FaceColor', 'k', 'EdgeColor', 'k'); ylabel('# CV folds'); xlabel('Cross validated Deviance for GLM (leave one vocalization out)');
        hold on; MeanLine = vline(mean(ModelGLMAmpDev_cv{nc})); MeanLine.LineWidth=2; MeanLine.LineStyle = '-';
        hold on;STELine1=vline(mean(ModelGLMAmpDev_cv{nc})-2*std(ModelGLMAmpDev_cv{nc})/NStims^0.5); STELine1.LineWidth=2;
        hold on;STELine2=vline(mean(ModelGLMAmpDev_cv{nc})+2*std(ModelGLMAmpDev_cv{nc})/NStims^0.5);STELine2.LineWidth=2;
        
        figure(6)
         subplot(2,2,3);hold on; histogram(ModelCeilR2_cv{nc}, 'FaceColor', [0 0.447 0.741], 'EdgeColor', [0 0.447 0.741]); 
        hold on; MeanLine = vline(mean(ModelCeilR2_cv{nc})); MeanLine.LineWidth=2; MeanLine.LineStyle = '-';MeanLine.Color = 'k';
        hold on;STELine1=vline(mean(ModelCeilR2_cv{nc})-2*std(ModelCeilR2_cv{nc})/NStims^0.5); STELine1.LineWidth=2; STELine1.Color = 'k';
        hold on;STELine2=vline(mean(ModelCeilR2_cv{nc})+2*std(ModelCeilR2_cv{nc})/NStims^0.5);STELine2.LineWidth=2;STELine2.Color = 'k';
        subplot(2,2,4);hold on; histogram(ModelGLMCeilDev_cv{nc}, 'FaceColor', [0 0.447 0.741], 'EdgeColor', [0 0.447 0.741]); 
        hold on; MeanLine = vline(mean(ModelGLMCeilDev_cv{nc})); MeanLine.LineWidth=2; MeanLine.LineStyle = '-';MeanLine.Color = 'k';
        hold on;STELine1=vline(mean(ModelGLMCeilDev_cv{nc})-2*std(ModelGLMCeilDev_cv{nc})/NStims^0.5); STELine1.LineWidth=2; STELine1.Color = 'k';
        hold on;STELine2=vline(mean(ModelGLMCeilDev_cv{nc})+2*std(ModelGLMCeilDev_cv{nc})/NStims^0.5);STELine2.LineWidth=2;STELine2.Color = 'k';
        
        subplot(2,2,4);hold on; histogram(ModelGLMFloorDev_cv{nc}, 'FaceColor', [0.494 0.184 0.556], 'EdgeColor', [0.494 0.184 0.556]); 
        hold on; MeanLine = vline(mean(ModelGLMFloorDev_cv{nc})); MeanLine.LineWidth=2; MeanLine.LineStyle = '-';MeanLine.Color = 'k';
        hold on;STELine1=vline(mean(ModelGLMFloorDev_cv{nc})-2*std(ModelGLMFloorDev_cv{nc})/NStims^0.5); STELine1.LineWidth=2; STELine1.Color = 'k';
        hold on;STELine2=vline(mean(ModelGLMFloorDev_cv{nc})+2*std(ModelGLMFloorDev_cv{nc})/NStims^0.5);STELine2.LineWidth=2;STELine2.Color = 'k';
          
        ModelGLMAmpGoodnessOfFit(nc,1) = mean((ModelGLMFloorDev_cv{nc} - ModelGLMAmpDev_cv{nc}) ./ (ModelGLMFloorDev_cv{nc} - ModelGLMCeilDev_cv{nc}));
        ModelGLMAmpGoodnessOfFit(nc,2) = std((ModelGLMFloorDev_cv{nc} - ModelGLMAmpDev_cv{nc}) ./ (ModelGLMFloorDev_cv{nc} - ModelGLMCeilDev_cv{nc}))/(length(ModelGLMCeilDev_cv{nc}))^0.5;
        
        % Full model to get the optimal filter
        ModelAmp = fitlm(AmpScore(:,1:NPC), Y_ZS );
        [ModelGLMAmpB, ModelGLMAmpDev, ModelGLMAmpStats] = glmfit(AmpScore(:,1:NPC), Y, 'Poisson');
        if ModelAmp.Rsquared.Adjusted>0
            ModelAmpLM_Coeff{nc,2} =  [ModelAmp.Coefficients.Estimate(1)' [ModelAmp.Coefficients.Estimate(2:(NPC+1))' zeros(1,size(XAmp,2)-NPC)]* AmpPC];
            ModelAmpLM_Coeff{nc,1} = -(Bins*TR - BestShift(cc));
        end
        ModelAmpGLM_Coeff{nc,2} =  [ModelGLMAmpB(1)' [ModelGLMAmpB(2:(NPC+1))' zeros(1,size(XAmp,2)-NPC)]* AmpPC];
        ModelAmpGLM_Coeff{nc,1} = -(Bins*TR - BestShift(cc));
        ModelAmpR2(nc,3) = ModelAmp.Rsquared.Adjusted;
        AnovaAmp = anova(ModelAmp, 'summary');
        ModelAmpR2(nc,1) = AnovaAmp.F(contains(AnovaAmp.Properties.RowNames, 'Model'));
        ModelAmpR2(nc,2) = AnovaAmp.pValue(contains(AnovaAmp.Properties.RowNames, 'Model'));
        figure(6);subplot(2,2,1); plot(Cell.MotorCoherenceOperant.CoherencyT_xTimeDelay, Cell.MotorCoherenceOperant.CoherencyT_filt, 'LineWidth',2);xlabel('Time Delay (ms)');ylabel('Coherency')
        Xlim = get(gca, 'XLim');
        subplot(2,2,2); plot(ModelAmpLM_Coeff{nc,1}, ModelAmpLM_Coeff{nc,2}(2:end), 'LineWidth',2); xlabel('Time Delay (X-Y) (ms)'); ylabel('Amplitude Linear Model coefficients (Y~Xamp)');
        hold on; yyaxis right;plot(ModelAmpGLM_Coeff{nc,1},ModelAmpGLM_Coeff{nc,2}(2:end), 'LineWidth',2);ylabel('Amplitude Poisson Generalized Linear Model coefficients (Y~Xamp)');
        set(gca, 'XLim', Xlim)
        subplot(2,2,3)
        title(sprintf('Adjusted R2 = %.3f CVR2 = %.3f pVal of F test vs constant model= %.2f', ModelAmp.Rsquared.Adjusted,mean(ModelAmpR2_cv{nc}), AnovaAmp.pValue(contains(AnovaAmp.Properties.RowNames, 'Model'))));
        subplot(2,2,4)
        title(sprintf('GLM goodness of fit = %.2f +/- %.2f', ModelGLMAmpGoodnessOfFit(nc,1), ModelGLMAmpGoodnessOfFit(nc,2)))
         % The fit is constrained by the variability in the neural response that we cannot estimate based on the data (no exact repetition of the same stimulus) however
         % we can model the neural noise as that of a Poisson process and
         % calculate what portion of the variance is left to explained by acoustic features -> what is the actual best R2 we could ever get given the data.
         % Bootstrap the calculation of random spike patterns that follow a
         % time varying Poisson of Lambda = Y, and calculate for each the
         % R2 (1- SSEm/SSE0)
         % This is the null model sum of squared errors
         SSE0 = sum((Y_ZS - mean(Y_ZS)).^2);
         % Generate random Poisson spike trains that follow the time
         % varying rate in YPerStim
         [Spike_array, Spike_arrival] = poisson_spiketrain_gen(YPerStim, YPerStimt,TLim, TR, Response_samprate, NBoots);
         
         
         % Run through the number of bootstrap loops to calculate ceil R2 and ceil deviance
         LL_Saturatedbb = LL_Calculus(YCount,YCount);
         for tt=1:size(Spike_arrival,2)
             [YPerStimbb, YPerStimtbb,~,YCountPerStimbb] = get_y_4Coherence(Spike_arrival(:,tt), Cell.Duration(IndVoc),TLim,TR);
             if InputFig
                 for Stim = 1:NStims
                     figure(16);
                     clf
                     Win = YPerStimt{Stim} + TR - TLim(1);
                     hist(Spike_arrival{Stim,tt},Win(end))
                     title(sprintf('NSpike = %d, SSE0 = %.2f SSECeil = %.2f',length(Spike_arrival{Stim,tt}), sum((YPerStim{Stim}-mean(Y)).^2), sum((YPerStim{Stim}-YPerStimbb{Stim}).^2)))
                     xlim([0 Win(end)]+TLim(1))
                     h=findobj(gca, 'Type', 'patch');
                     h.FaceColor = 'k';
                     h.EdgeColor = 'k';
                     hold on
                     yyaxis right
                     plot(YPerStimt{Stim},YPerStim{Stim}, 'LineWidth',2)
                     ylabel('Original rate mHz')
                     hold on
                     plot(YPerStimtbb{Stim},YPerStimbb{Stim}, 'LineWidth',2)
                     hold on
                     Hmean=hline(mean(Y));
                     Hmean.LineWidth = 2;
                     pause
                 end
             end    
             Y_bb = ([YPerStimbb{:}]' - Y_Mu)/Y_Sigma;
             SSECeil(nc,tt) = sum((Y_ZS - Y_bb).^2);
             R2CeilBoot(nc,tt) = 1-SSECeil(nc,tt)/SSE0;
             YCountbb = [YCountPerStimbb{:}]';
             LL_Ceilbb = LL_Calculus(YCount, YCountbb);
             DevianceCeilBoot(nc,tt) = -2*(LL_Ceilbb - LL_Saturatedbb);
         end
         figure(6)
         subplot(2,2,3);hold on; histogram(R2CeilBoot(nc,:), 'FaceColor', 'g', 'EdgeColor', 'g'); 
        hold on; MeanLine = vline(mean(R2CeilBoot(nc,:))); MeanLine.LineWidth=2; MeanLine.LineStyle = '-';MeanLine.Color = 'k';
        hold on;STELine1=vline(mean(R2CeilBoot(nc,:))-2*std(R2CeilBoot(nc,:))/NStims^0.5); STELine1.LineWidth=2; STELine1.Color = 'k';
        hold on;STELine2=vline(mean(R2CeilBoot(nc,:))+2*std(R2CeilBoot(nc,:))/NStims^0.5);STELine2.LineWidth=2;STELine2.Color = 'k';
        subplot(2,2,4);hold on; histogram(DevianceCeilBoot(nc,:), 'FaceColor', 'g', 'EdgeColor', 'g'); 
        hold on; MeanLine = vline(mean(DevianceCeilBoot(nc,:))); MeanLine.LineWidth=2; MeanLine.LineStyle = '-';MeanLine.Color = 'k';
        hold on;STELine1=vline(mean(DevianceCeilBoot(nc,:))-2*std(DevianceCeilBoot(nc,:))/NStims^0.5); STELine1.LineWidth=2; STELine1.Color = 'k';
        hold on;STELine2=vline(mean(DevianceCeilBoot(nc,:))+2*std(DevianceCeilBoot(nc,:))/NStims^0.5);STELine2.LineWidth=2;STELine2.Color = 'k';
             
        
        % revise this does not work with actual model with several time
        % points
        XPredictAmp = min(XAmp):(max(XAmp)-min(XAmp))/10:max(XAmp);
        XPredictAmp_ZS = (XPredictAmp-XAmp_Mu)./XAmp_Sigma;
        [YPredictAmp_ZS,YPredictAmpci_ZS] = predict(ModelAmp,XPredictAmp_ZS');
        
        ModelwoSal = fitlm(TableFull, 'Rate ~ Amplitude + SpectralMean + CallType +  Amplitude:SpectralMean + Amplitude:CallType +  SpectralMean:CallType' , 'CategoricalVars', 4);
        ModelSal = fitlm(XSal,Y_ZS);
        ModelSalR2(nc,3) = ModelSal.Rsquared.Adjusted;
        AnovaSal = anova(ModelSal, 'summary');
        ModelSalR2(nc,1) = AnovaSal.F(contains(AnovaSal.Properties.RowNames, 'Model'));
        ModelSalR2(nc,2) = AnovaSal.pValue(contains(AnovaSal.Properties.RowNames, 'Model'));
        if ModelSal.Rsquared.Adjusted>0
            ModelSalR2(nc,4) =  ModelSal.Coefficients.Estimate(2);
        end
        XPredictSal = min(XSal):(max(XSal)-min(XSal))/10:max(XSal);
        XPredictSal_ZS = (XPredictSal-XSal_Mu)./XSal_Sigma;
        [YPredictSal_ZS,YPredictSalci_ZS] = predict(ModelSal,XPredictSal_ZS');
        
        ModelwoSpecMean = fitlm(TableFull, 'Rate ~ Amplitude + Saliency + CallType + Amplitude:Saliency +  Amplitude:CallType +  Saliency:CallType' , 'CategoricalVars', 4);
        ModelSpecMean = fitlm(XSpecMean,Y_ZS);
        ModelSpecMeanR2(nc,3) = ModelSpecMean.Rsquared.Adjusted;
        AnovaSpecMean = anova(ModelSpecMean, 'summary');
        ModelSpecMeanR2(nc,1) = AnovaSpecMean.F(contains(AnovaSpecMean.Properties.RowNames, 'Model'));
        ModelSpecMeanR2(nc,2) = AnovaSpecMean.pValue(contains(AnovaSpecMean.Properties.RowNames, 'Model'));
        if ModelSpecMean.Rsquared.Adjusted>0
            ModelSpecMeanR2(nc,4) =  ModelSpecMean.Coefficients.Estimate(2);
        end
        XPredictSpecMean = min(XSpecMean):(max(XSpecMean)-min(XSpecMean))/10:max(XSpecMean);
        XPredictSpecMean_ZS = (XPredictSpecMean-XSpecMean_Mu)./XSpecMean_Sigma;
        [YPredictSpecMean_ZS,YPredictSpecMeanci_ZS] = predict(ModelSpecMean,XPredictSpecMean_ZS');
        
        ModelwoCallType = fitlm(TableFull, 'Rate ~ Amplitude + Saliency + SpectralMean +  Amplitude:Saliency + Amplitude:SpectralMean + Saliency:SpectralMean ' , 'CategoricalVars', 4);
        ModelCallType = fitlm(categorical(Trill1_Bark0_local),Y_ZS);
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
        
        figure(4)
        clf
        subplot(1,3,1)
        shadedErrorBar(XPredictAmp,(YPredictAmp_ZS.*Y_Sigma + Y_Mu)/(10^-3),([YPredictAmp_ZS-YPredictAmpci_ZS(:,1) YPredictAmpci_ZS(:,2)- YPredictAmp_ZS]'.*Y_Sigma + Y_Mu)./(10^-3))
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
        shadedErrorBar(XPredictSal,(YPredictSal_ZS.*Y_Sigma + Y_Mu)/(10^-3),([YPredictSal_ZS-YPredictSalci_ZS(:,1) YPredictSalci_ZS(:,2)- YPredictSal_ZS]'.*Y_Sigma + Y_Mu)./(10^-3))
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
        shadedErrorBar(XPredictSpecMean,(YPredictSpecMean_ZS.*Y_Sigma + Y_Mu)/(10^-3),([YPredictSpecMean_ZS-YPredictSpecMeanci_ZS(:,1) YPredictSpecMeanci_ZS(:,2)- YPredictSpecMean_ZS]'.*Y_Sigma + Y_Mu)./(10^-3))
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
        
        suplabel(sprintf('Cell %d/%d   %s   TR = %d ms   R2 = %.2f   Info = %.2f bits/s    pRandSpike = %.2f', nc, NCells,MC.CellsPath(cc).name(1:end-4), TR,FullModelR2(nc,3), Cell.MotorCoherenceOperant.Info, Cell.MotorCoherenceOperant.Info_pRandSpikePerm), 't');
        suplabel(sprintf('CallType R2 = %.2f Partial R2 = %.2f p=%.2f', ModelCallTypeR2(nc,3), CallTypePredictor(nc,[3 2])),'x');
        
    
    pause(1)
    
end
CellsPath = MC.CellsPath;
save(fullfile(HDPath, 'LM_Acoustic.mat'),'CellsPath','GoodInfo','NCells', 'AmpPredictor','SalPredictor','SpecMeanPredictor' , 'CallTypePredictor','FullModelR2','ModelAmpR2','ModelSalR2','ModelSpecMeanR2','ModelCallTypeR2');
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
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,sprintf('SU =%d',sum(SU1MU0(GoodInfo))),'Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,sprintf('MU =%d',sum(~SU1MU0(GoodInfo))),'Color',[0 0 0], 'FontWeight','bold')

subplot(2,2,2)
scatter(Info(GoodInfo), FullModelR2(:,3), 30, [FullModelR2(:,2)<0.01 zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlabel('Info on Coherence with Amplitude (bits)')
ylabel('Adjusted R2 full linear Model')
xlim([0 5])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,sprintf('Significant (N=%d)', sum(FullModelR2(:,2)<0.01)),'Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,sprintf('Non-Sig (N=%d)', sum(FullModelR2(:,2)>=0.01)),'Color',[0 0 0], 'FontWeight','bold')

subplot(2,2,3)
histogram(FullModelR2(:,3), 'BinWidth',0.01,'FaceColor','k')
ylabel('# High Info Cells')
xlabel('Adjusted R2 full linear model')
suplabel(sprintf('Full acoustic model performance N=%d cells', length(GoodInfo)),'t')

% Now figure of individual models with all cells
figure(32)
clf
subplot(3,4,1)
scatter(FullModelR2(:,3), ModelAmpR2(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hold on
plot([-0.2 0.4], [-0.2 0.4],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Amplitude')

subplot(3,4,2)
scatter(FullModelR2(:,3), ModelSalR2(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hold on
plot([-0.2 0.4], [-0.2 0.4],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Pitch Saliency')

subplot(3,4,3)
scatter(FullModelR2(:,3), ModelSpecMeanR2(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hold on
plot([-0.2 0.4], [-0.2 0.4],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 Spectral Mean')

subplot(3,4,4)
scatter(FullModelR2(:,3), ModelCallTypeR2(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hold on
plot([-0.2 0.4], [-0.2 0.4],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('R2 Full Model')
ylabel('R2 CallType')

subplot(3,4,5)
scatter(FullModelR2(:,3), AmpPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hold on
plot([-0.2 0.4], [-0.2 0.4],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Amplitude')

subplot(3,4,6)
scatter(FullModelR2(:,3), SalPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hold on
plot([-0.2 0.4], [-0.2 0.4],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,7)
scatter(FullModelR2(:,3), SpecMeanPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hold on
plot([-0.2 0.4], [-0.2 0.4],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Full Model R2')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,8)
scatter(FullModelR2(:,3), CallTypePredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hold on
plot([-0.2 0.4], [-0.2 0.4],'k--', 'LineWidth',2)
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
ylabel('Partial R2 Call Type')
xlabel('Full Model R2')

subplot(3,4,9)
scatter(ModelAmpR2(:,4), AmpPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 0.4])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
hold on
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Amplitude')
ylabel('Partial R2 Amplitude')

subplot(3,4,10)
scatter(ModelSalR2(:,4), SalPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 0.4])
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
hold on

text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Pitch Saliency')
ylabel('Partial R2 Pitch Saliency')

subplot(3,4,11)
scatter(ModelSpecMeanR2(:,4), SpecMeanPredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 0.4])
hold on
YLim = get(gca, 'YLim');
XLim = get(gca, 'XLim');
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9.5/10,'SU','Color',[1 0 0], 'FontWeight','bold')
text(diff(XLim)*3/4 + XLim(1), YLim(2)*9/10,'MU','Color',[0 0 0], 'FontWeight','bold')
xlabel('Coefficient Spectral Mean')
ylabel('Partial R2 Spectral Mean')

subplot(3,4,12)
scatter(ModelCallTypeR2(:,4), CallTypePredictor(:,3),30,[SU1MU0(GoodInfo) zeros(size(GoodInfo)) zeros(size(GoodInfo))], 'filled')
ylim([-0.2 0.4])
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
    DefaultValue = NaN;
else
    DefaultValue = DefaultVal;
end

for stim = 1:length(Duration)
    if abs(round(length(BioSound{stim}.sound)/BioSound{stim}.samprate*10^3) - round(Duration(stim)))>1
        keyboard
    end
    % Get ready an output vector for the stim acoustic features that was sampled at 1000Hz
    Xt = Delay(1) : (Delay(2) + Duration(stim));
    XPerStim_temp = DefaultValue.*ones(size(Xt));
    if any(Xt>=0)
        FeatureVal = BioSound{stim}.(sprintf('%s',Feature));
        Xtpos = find(Xt>=0);
        XPerStim_temp(Xtpos(1:min(length(Xtpos), length(FeatureVal)))) = FeatureVal(1:min(length(Xtpos), length(FeatureVal)));
    end  
    
    if ~ischar(DefaultVal)
        XPerStim_temp(isnan(XPerStim_temp)) = DefaultValue;
    else
        XPerStim_temp(isnan(XPerStim_temp)) = nanmean(FeatureVal);
    end
    % Warning, resample gives very weird results here!!!
    %     [XPerStim{stim}] = resample(XPerStim_temp, Fs, 1000);
    %     [XPerStim{stim},XPerStimt{stim}] = resample(XPerStim_temp, (1:length(XPerStim_temp)),Fs*10^-3, 'spline');
    
    % I'm doing my own resampling
    % Time slots
    TimeBinsXOnsetInd = round(1 :round(TR-Overlap): (Delay(2) - Delay(1) + Duration(stim))); % These are slightly different than in get_Y_4GLM, because the times slot are used as indices in the vector and not as actuel time values!
    TimeBinsXOffsetInd = TimeBinsXOnsetInd + TR -1;
    TimeBinsXOnsetInd = TimeBinsXOnsetInd(TimeBinsXOffsetInd<=(Delay(2) - Delay(1) + Duration(stim))); % Only keep windows that are within the call
    TimeBinsXOffsetInd = TimeBinsXOffsetInd(TimeBinsXOffsetInd<=(Delay(2) - Delay(1) + Duration(stim))); % Only keep windows that are within the call
    
    
    TimeBinsXOnset = Delay(1) :round(TR-Overlap): (Delay(2) + Duration(stim));
    TimeBinsXOffset = TimeBinsXOnset + TR;
    TimeBinsXOnset = TimeBinsXOnset(TimeBinsXOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    TimeBinsXOffset = TimeBinsXOffset(TimeBinsXOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
    XPerStimt{stim} = TimeBinsXOnset + (TimeBinsXOffset - TimeBinsXOnset)/2;
    
    XPerStim_resamp = nan(TR, length(TimeBinsXOnsetInd));
    for tt=1:length(TimeBinsXOnsetInd)
        XPerStim_resamp(:,tt) = XPerStim_temp(TimeBinsXOnsetInd(tt): TimeBinsXOffsetInd(tt))';
    end
    XPerStim{stim} = mean(XPerStim_resamp);
    
    
    if DebugFig
        figure(200)
        clf
        plot(Xt,XPerStim_temp, 'LineWidth',2)
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

function [XPerStim] = get_x_4gammafit(BioSound, Duration, SATPerStim, FilterSize, DefaultVal,Feature)
% This function calculate for each stim and at each spike arrival time SAT,
% the value of amplitude in a window FilterSize centered around SAT.
DebugFig=0;

% calculate the cell array of acoustic data for each cell. For each
% voc, XPerStim is a matrix where each column is a spike and each row corresponds to the values of
% the acoustic feature in the window starting FilterSize/2 ms before SAT and FilterSize/2 ms after.
%  Fs = round(1/((TR-Overlap).*10^-3));
% Vectors of the acoustic features
XPerStim = cell(1,length(Duration));
if ischar(DefaultVal) && strcmp(DefaultVal, 'mean')
    DefaultValue = NaN;
else
    DefaultValue = DefaultVal;
end

for stim = 1:length(Duration)
    if abs(round(length(BioSound{stim}.sound)/BioSound{stim}.samprate*10^3) - round(Duration(stim)))>1
        keyboard
    end
    % Get the time varying feature for that stim 
    FeatureVal = BioSound{stim}.(sprintf('%s',Feature));
    % loop through spikes and get the corresponding feature value vector
    XPerStim{stim} = DefaultValue .* ones(FilterSize+1, length(SATPerStim{stim}));
    for isp = 1:length(SATPerStim{stim})
        SAT = SATPerStim{stim}(isp);
        Xt = round(SAT + (-FilterSize/2 : FilterSize/2));
        if any(Xt>=0)
            Xtpos = find(Xt>=0);
            XPerStim{stim}(Xtpos(1:min(length(Xtpos), length(FeatureVal))), isp) = FeatureVal(1:min(length(Xtpos), length(FeatureVal)));
        end  
    end
    
    
    if DebugFig
        figure(200)
        clf
        imagesc(XPerStim{stim}'); colorbar()
        set(gca, 'XTick', (1 :25: FilterSize+1),'XTickLabel', (-FilterSize/2 :25: FilterSize/2))
        xlabel('Time ms')
        ylabel('spike #')
        title(sprintf('%s Stim %d/%d',Feature, stim,length(Duration)));
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

function [YPerStim, YPerStimt,SpikePattern,YCountPerStim] = get_y_4Coherence(SAT, Duration,Delay,TR,Overlap, Fs)
DebugFig = 0;
if nargin<5
    Overlap = 0;
end
if length(Delay)==1
    Delay = [-Delay Delay];
end
% Calculate the time varying rate applying a gaussian window TR on the
% spike pattern. The spike pattern considered starts -Delay ms
% before the onset of the vocalization and stops Delay ms after the
% offset of the vocalization
YPerStim = cell(1,length(Duration));
YPerStimt = cell(1,length(Duration));
YCountPerStim = cell(1,length(Duration));
% Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
nStd =(max(Duration) - Delay(1) + Delay(2))/10; % before set as 4
Tau = (TR/2);
T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
Expwav = Expwav./sum(Expwav);
% Frequency at which the neural data should be sampled
if nargin<6
    Fs = round(1/((TR-Overlap).*10^-3));
end
SpikePattern = cell(length(Duration),1);
% Loop through the stimuli and fill in the matrix
for stim=1:length(Duration)
    
    % Time slots for the neural response
    TimeBinsY = Delay(1) : (Delay(2) + round(Duration(stim)));
    SpikePattern{stim} = zeros(1,length(TimeBinsY)-1);
    for isp = 1:length(SAT{stim})
        SpikeInd = round(SAT{stim}(isp));
        if (SpikeInd>=Delay(1)) && (SpikeInd<(Delay(2) + Duration(stim)))
            SpikePattern{stim}(SpikeInd - Delay(1) +1) = SpikePattern{stim}(SpikeInd - Delay(1) +1) +1;
        end
    end
    
    % Convolve with Gaussian to obtain our smooth time varying spike train
    % and resample if necessary
    if Fs == 1000
        YPerStim{stim} = conv(SpikePattern{stim}, Expwav,'same');
        YPerStimt{stim} = TimeBinsY;
    else
        YPerStim_local = conv(SpikePattern{stim}, Expwav,'same');
        if sum(YPerStim_local)>0
            YPerStim_local = YPerStim_local/sum(YPerStim_local)*sum(SpikePattern{stim}); % Make sure we keep the right number of sipkes after convolution!
        end
        
        % resampling function is really doing weird things at edges...
        % doing my own resampling
        TimeBinsYOnsetInd = round(1 :round(TR-Overlap): (Delay(2) - Delay(1) + Duration(stim))); % These are slightly different than in get_Y_4GLM, because the times slot are used as indices in the vector and not as actuel time values!
        TimeBinsYOffsetInd = TimeBinsYOnsetInd + TR -1;
        TimeBinsYOnsetInd = TimeBinsYOnsetInd(TimeBinsYOffsetInd<=(Delay(2) - Delay(1) + Duration(stim))); % Only keep windows that are within the call
        TimeBinsYOffsetInd = TimeBinsYOffsetInd(TimeBinsYOffsetInd<=(Delay(2) - Delay(1) + Duration(stim))); % Only keep windows that are within the call
        
        
        TimeBinsYOnset = Delay(1) :round(TR-Overlap): (Delay(2) + Duration(stim));
        TimeBinsYOffset = TimeBinsYOnset + TR;
        TimeBinsYOnset = TimeBinsYOnset(TimeBinsYOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
        TimeBinsYOffset = TimeBinsYOffset(TimeBinsYOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
        YPerStimt{stim} = TimeBinsYOnset + (TimeBinsYOffset - TimeBinsYOnset)/2; 
        
        YPerStim_resamp = nan(TR, length(TimeBinsYOnsetInd));
        YCountPerStim{stim} = nan(1,length(TimeBinsYOnsetInd));
        for tt=1:length(TimeBinsYOnsetInd)
            YPerStim_resamp(:,tt) = YPerStim_local(TimeBinsYOnsetInd(tt): TimeBinsYOffsetInd(tt))';
            YCountPerStim{stim}(tt) = sum(SpikePattern{stim}(TimeBinsYOnsetInd(tt): TimeBinsYOffsetInd(tt)));
        end
        YPerStim{stim} = mean(YPerStim_resamp);
        
        
        
        if DebugFig
            figure(200)
            clf
            plot((Delay(1)+0.5):(Duration(stim)+Delay(2)),YPerStim_local, 'LineWidth',2)
            xlabel('Time ms')
            ylabel('Spike Rate mHz (/ms)')
            title(sprintf('Stim %d/%d',stim,length(Duration)));
            %         RemainTime = length(XPerStim_temp) - (length(XPerStim{stim})-1)*(1/Fs*10^3);
            %         XPerStimt{stim} = RemainTime/2+(1/Fs*10^3)*(0:(length(XPerStim{stim})-1));
            hold on
            plot(YPerStimt{stim},YPerStim{stim}, 'LineWidth',2)
            legend({'original' 'resampled'}, 'AutoUpdate','off')
            hold on
            SpikeTimes = TimeBinsY(logical(SpikePattern{stim}));
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



function [YPerStim, YPerStimt,SpikePattern, SATPerStim, ISIPerStim] = get_y_4GammaFit(SAT, Duration,Delay,TR, ONNK)
if nargin<5
    ONNK = 4;
end
DebugFig=0;
if length(Delay)==1
    Delay = [-Delay Delay];
end
% Calculate the time varying rate applying a time varying gaussian window TR on the
% spike pattern. The spike pattern considered starts -Delay ms
% before the onset of the vocalization and stops Delay ms after the
% offset of the vocalization
YPerStim = cell(1,length(Duration));
YPerStimt = cell(1,length(Duration));
SATPerStim = cell(1,length(Duration));
ISIPerStim = cell(1,length(Duration));
% Gaussian window of std equal to Tau (68% of Gaussian centered in 2*Tau)
% Tau is chosen as the 4th closest neighbor spike
nStd =(max(Duration) - Delay(1) + Delay(2))/10; % before set as 4


SpikePattern = cell(length(Duration),1);
% Loop through the stimuli and fill in the matrix
for stim=1:length(Duration)
    
    % Time slots for the neural response
    TimeBinsY = Delay(1) : (Delay(2) + round(Duration(stim)));
    % Spike pattern for that stim
    SpikePattern{stim} = zeros(1,length(TimeBinsY)-1);
    % make sure Spike arrival times are chronological
    SAT_local = sort(SAT{stim});
    % get the indices of spikes within that stim
    SpikeInd = round(SAT_local);
    if length(SpikeInd)<(ONNK+1)
        warning('Too few spike for this event no calculation')
        continue
    end
    SpikeInd = find((SpikeInd>=Delay(1)) .* (SpikeInd<(Delay(2) + Duration(stim))));
    NSpikes = length(SpikeInd);
    SATPerStim{stim} = SAT_local(SpikeInd);
    ISIPerStim{stim} = nan(1,length(SpikeInd));
    % Now for each spike fillin the spike pattern and convolve it in NNKE with a gaussian which std is taken to
    % be the 4th nearest neighbor (ONNK)
    NNKE = zeros(NSpikes,length(TimeBinsY));
    GoodSpikes = ones(NSpikes,1);
    for isp = 1:NSpikes
        SpikePattern{stim}(round(SATPerStim{stim}(isp)) - Delay(1) +1) = SpikePattern{stim}(round(SATPerStim{stim}(isp)) - Delay(1) +1) +1;
        NNKE(isp,(round(SATPerStim{stim}(isp)) - Delay(1) +1)) = 1;
        % Now try a time varying kernel for rate estimation (nearest neighbor?)
        if (isp-ONNK)>0 && (isp+ONNK)<=length(SATPerStim{stim})
            Tau = min(abs(SATPerStim{stim}(isp) - [SATPerStim{stim}(isp+ONNK) SATPerStim{stim}(isp-ONNK)]));
        elseif (isp-ONNK)>0 && (SpikeInd(isp) + ONNK)<=length(SAT_local) % The ONNKthnearest spike in the future is out of the vocalization zone but still available
            Tau = min(abs(SATPerStim{stim}(isp)- [SATPerStim{stim}(isp-ONNK) SAT_local(SpikeInd(isp) + ONNK)]) );
        elseif (isp+ONNK)<=length(SATPerStim{stim}) && (SpikeInd(isp)-ONNK)>0 % The ONKth nearest spike in the past is out of the vocalization zone but still available
            Tau = min(abs(SATPerStim{stim}(isp)- [SATPerStim{stim}(isp+ONNK) SAT_local(SpikeInd(isp) - ONNK)]));
        elseif (isp-ONNK)>0 % There is no accessible ONNKth nearest spike in the future only take the distance to the ONNKth spike in the past
            Tau = abs(SATPerStim{stim}(isp)- SATPerStim{stim}(isp-ONNK));
        elseif (isp+ONNK)<=length(SATPerStim{stim}) % there is no accessible ONNKth nearest spike in the past, only take distance to the OKKNth spike in the future
            Tau = abs(SATPerStim{stim}(isp)- SATPerStim{stim}(isp+ONNK));
        elseif (SpikeInd(isp)-ONNK)>0 && (SpikeInd(isp) + ONNK)<=length(SAT_local) % There are accessible spikes outside of the vocalization zone both in future and past
            Tau = min(abs(SATPerStim{stim}(isp)- [SAT_local(SpikeInd(isp) - ONNK) SAT_local(SpikeInd(isp) + ONNK)]) );
        elseif (SpikeInd(isp)-ONNK)>0  % There are accessible spikes outside of the vocalization zone only in past
            Tau = abs(SATPerStim{stim}(isp)- SAT_local(SpikeInd(isp) - ONNK));
        elseif (SpikeInd(isp) + ONNK)<=length(SAT_local)  % There are accessible spikes outside of the vocalization zone only in future
            Tau = abs(SATPerStim{stim}(isp)- SAT_local(SpikeInd(isp) + ONNK));
        elseif length(SAT_local)<(2*ONNK+1) % this is expected take a default value...
            Tau = max(Tau,TR/2);
        else
            warning('ISSUE CANNOT FIND THE %dth nearest spike!!', ONNK);
            keyboard
        end
        Tau = max(Tau,TR/2);
        
        % Now get the time to the previous spike for each spike
        if (isp-1)>0
            ISIPerStim{stim}(isp) = abs(SATPerStim{stim}(isp) - SATPerStim{stim}(isp-1));
        elseif (SpikeInd(isp)-1)>0 % The 1st nearest spike in the past is out of the vocalization zone but still available
            ISIPerStim{stim}(isp)  = abs(SATPerStim{stim}(isp)- SAT_local(SpikeInd(isp) - 1));
        else
            warning('Cannot find a previous spike cancelling that spike')
            GoodSpikes(isp) = 0;
        end
        % construct the Gaussian
        T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
        Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
        Expwav = Expwav./sum(Expwav);
        NNKE(isp,:) = conv(NNKE(isp,:), Expwav,'same');
    end
    ISIPerStim{stim} = ISIPerStim{stim}(logical(GoodSpikes));
    SATPerStim{stim} = SATPerStim{stim}(logical(GoodSpikes));
    NNKE = NNKE(logical(GoodSpikes),:);
    YPerStim{stim} = sum(NNKE,1);
    if length(YPerStim{stim})~= length(TimeBinsY)
        warning('Unexpected length')
        keyboard
    end
    YPerStimt{stim} = TimeBinsY;
    if DebugFig
        figure(200); clf
        hist(SATPerStim{stim}); hold on; yyaxis right; plot(YPerStimt{stim}, YPerStim{stim}, '-r','LineWidth',2); xlabel('Time (ms)')
        title(sprintf('%dth nearest neighbor rate estimation STim %d/%d', ONNK, stim,length(Duration)));
        pause(1)
    end
    
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
    
    xlim([min(Delay(1),YPerStimt{stim}(1)) max(Duration(stim) + Delay(2), YPerStimt{stim}(end))])
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
    line([YPerStimt{stim}(1) YPerStimt{stim}(end)], [0.875 0.875].*YLIM(2), 'Color',ColorCode(1,:), 'LineWidth',12)
    hold on
    text((YPerStimt{stim}(end)-YPerStimt{stim}(1))/2.3 +YPerStimt{stim}(1) ,0.875.*YLIM(2),'Y', 'Color', [1 1 1])
    hold on
    line([XPerStimt{stim}(1) XPerStimt{stim}(end)], [0.775 0.775].*YLIM(2), 'Color',ColorCode(2,:), 'LineWidth',12)
    hold on
    text((XPerStimt{stim}(end)-XPerStimt{stim}(1))/2.3 + XPerStimt{stim}(1),0.775.*YLIM(2),'X', 'Color', [1 1 1])
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
        ylim([0 YLim(2)*1.5])
    end
    xlim([min(Delay(1),YPerStimt{stim}(1)) max(Duration(stim) + Delay(2), YPerStimt{stim}(end))])
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


function [Spike_array, Spike_arrival] = poisson_spiketrain_gen(YPerStim, YPerStimt,TLim, TR, Response_samprate, NBoots)
if nargin<6
    NBoots = 100;
end
if nargin<5
    Response_samprate = 1000;
end
DebugFig=0;
% We generate spike trains (100 per stim) with inhomogenous Poisson distribution
NStims = length(YPerStim);
Spike_array = cell(NStims,NBoots); % contains the array of spikes for each stim (number of spikes in each 1ms bin for each trial)
YPerStim; % Ground truth Poisson rate: contains the time varying spike rate (# spikes /ms) for each stim that is used as an input mean to the inhomogenous Poisson responses
NSpikes = nan(NStims,1); % contains the total sum number of spikes for each stimulus accross the NTrials trials (one scalar per stimulus)
Spike_arrival = cell(NStims,NBoots); % contains all spike arrival times from all trials for each stim (one vector per stimulus)


for Stim=1:NStims
    Win = YPerStimt{Stim} + TR - TLim(1);% we want the end for each time bin and the time should be in reference to the onset of the Y estimation (-TLim)
    NSpikeslocal = nan(NBoots,1);
    Nb_Win = length(Win);
    
    for tt=1:NBoots
        Spike_array{Stim,tt}= zeros(Win(end)*Response_samprate/1000,1);
        NbSpikes_temp = poissrnd(YPerStim{Stim}.*TR);
        NSpikeslocal(tt) = sum(NbSpikes_temp);
        Win_old=0;
        for ww=1:(Nb_Win)
            if NbSpikes_temp(ww)>0
                Ind=[];
                NS_temp = NbSpikes_temp(ww);
                while NS_temp>0
                    Ind = [Ind randperm(TR*Response_samprate/1000,min(TR*Response_samprate/1000,NS_temp))];
                    NS_temp = NS_temp - min(TR*Response_samprate/1000,NS_temp);
                end
                for ii=1:NbSpikes_temp(ww)
                    Spike_array{Stim,tt}(Win_old*Response_samprate/1000+Ind(ii)) = Spike_array{Stim,tt}(Win_old*Response_samprate/1000+Ind(ii)) + 1;
                end
            end
            Win_old=Win(ww);
        end
    end
    
    % Total number of spikes over all trials of that stim
    NSpikes(Stim) =  sum(NSpikeslocal);
    
    
    % Getting spike arrival times of each trial for that stim
    
    
    for tt=1:NBoots
        Spike_arrival{Stim,tt} = nan(NSpikeslocal(tt),1);
        % Find the spike arival times for that trial and fill in
        % the by-trial set of spike arrival time. Refernce is stim
        % onset here
        TA = find(Spike_array{Stim,tt});
        spike_count = 0;
        for sp=1:length(TA)
            Ns_local = Spike_array{Stim,tt}(TA(sp));
            while Ns_local>0
                spike_count = spike_count +1;
                Spike_arrival{Stim,tt}(spike_count) = TA(sp)*1000/Response_samprate + TLim(1);
                Ns_local = Ns_local - 1;
            end
        end
    end
    
    
    % Plot the histogram of spike arrival times for the first botstrap for that stim
    if DebugFig
        figure(15);
        clf
        for tt=10:10:100
            subplot(10,1,tt/10)
            hist(Spike_arrival{Stim,tt},Win(end))
            title(sprintf('NSpike = %d',length(Spike_arrival{Stim,tt})))
            xlim([0 Win(end)]+TLim(1))
            h=findobj(gca, 'Type', 'patch');
            h.FaceColor = 'k';
            h.EdgeColor = 'k';
            hold on
            yyaxis right
            plot(YPerStimt{Stim},YPerStim{Stim}, 'LineWidth',2)
            ylabel('Original rate mHz')
        end
    end
    
end
end
