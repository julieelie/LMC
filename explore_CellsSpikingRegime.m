addpath(genpath('/Users/elie/Documents/CODE/operant_bats'))
addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'))
addpath(genpath('/Users/elie/Documents/CODE/LMC'))
addpath(genpath('/Users/elie/Documents/CODE/Kilosort2'))
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'))
addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))
Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';
% BasePath = '/Volumes/Julie4T';
BasePath = '/Volumes/server_home/users/JulieE/LMC/';

%% Generate the list of paths to gather the raw data
% OutputPath = '/Users/elie/Documents/LMCResults';
OutputPath = '/Volumes/server_home/users/JulieE/LMC/LMCResults';
[ListSSU] = gather_neural_datapath(BasePath);

%% Or alternatively use the already saved list of raw cells with the indices of cells that passed the sanitary check (= good single units) (done in wrapper_allbehav_neuro_script.m)
load('GoodCellIndicesAll.mat','GoodCellIndices','ListSSU')

%% Listing the cells for which the spiking activity has been alligned and extracted for vocal production and other behavioral events
AllFiles = dir(fullfile(OutputPath,'*.mat'));
Files2run = zeros(length(AllFiles),1);
for ff=1:length(AllFiles)
    if length(strfind(AllFiles(ff).name, '_'))==3
        Files2run(ff) = 1;
    end
end
CellsPath = AllFiles(logical(Files2run));

%% Explore the regime of one particularly interesting cell
NCells = length(GoodCellIndices);
NeuralBuffer = 5000; %used to extract neural data by neuralData_compile_per_file when called by wrapper_allbehav_neuro_script
NeuralWin = 2000; %duration of the time window in ms for the time varying analysis of the Regime
Overlap = 1000; %ms (step for sliding the NeuralWin
Buffer = 500; % if a vocalization is within Buffer ms of the edges of the NeuralWin window, then it is counted as during a vocalization
ss=510; %ss=211 (59834_20190617_TT3_SSS_155) %510 (59834_20190628_TT3_SSS_178) %450(59834_20190626_TT3_SSS_173) 210(59834_20190617_TT2_SSS_94) Counter example: 449 (T1)  (poisson everywhere? on T1) 212 and 451(gamma everywhere? both on T3)

% Have a look at the raster plot alligned to vocalization onset for that
% cell by running the following code (you want to stop the code plot_rastervoc_perfile at line 62
% to save the raster
% InputDataFile = ListSSU{GoodCellIndices(ss)};
% plot_rastervoc_perfile(InputDataFile, OutputPath, [2000 2000], 0, 0)

% Load the data
[Path2Cell,CellName] = fileparts(ListSSU{GoodCellIndices(ss)});
fprintf(1,'File %d/%d: %s \n',ss,NCells,CellName)
CellName_list{ss} = CellName;
Path2Cell = [BasePath Path2Cell(strfind(Path2Cell, 'LMC'):end)];
CellName_full = fullfile(Path2Cell, [CellName '.mat']);
load(CellName_full, 'Spike_arrival_times'); % load spike arrival times in useconds

Ind_ = strfind(CellName, '_');
SUBJ = CellName(1:(Ind_(1)-1));
Date = CellName(Ind_(1)+ (1:8));
Tetrode = CellName(Ind_(2)+ 3);
SSQ = CellName(Ind_(4) - 1);
CellID = CellName((Ind_(4) + 1) : end);
FileID = fullfile(OutputPath, [sprintf('%s_%s_SS%s_%s-%s',SUBJ, Date, SSQ, Tetrode, CellID) '.mat']);
load(FileID,'FreeBehavSession','OperantSession','PlayBackSession','QualitySSU','What','Who','ExpType','VocRank','Duration', 'DelayBefore','DelayAfter', 'SpikesArrivalTimes_Behav', 'SpikesArrivalTimes_Baseline') % these spike times are centered around vocalization onset


% Find vocalizations onset/offsets
PathParts = regexp(Path2Cell, '/', 'split');
Loggers_dir = ['/' fullfile(PathParts{1:end-2})];
VocDataFiles= dir(fullfile(Loggers_dir, sprintf('%s_*_VocExtractData.mat', Date(3:end))));
VocDataTime = cell(1,length(VocDataFiles));
for ff=1:length(VocDataFiles)
    load(fullfile(VocDataFiles(ff).folder, VocDataFiles(ff).name), 'Voc_transc_time_refined'); % this variable is in ms
    % Make sure to re-order values in chronological order
    [~,SortInd] = sort(Voc_transc_time_refined(:,1));
    VocDataTime{ff} = Voc_transc_time_refined(SortInd,:)';
end
VocDataTime = [VocDataTime{:}]'- OperantSession(1)*10^3; % Operant Session value is in s and VocDataTime in ms

% first extract the spike regime we can see on the raster plot for
% trill and non trill
TimeOn = -NeuralBuffer;
TimeOff = Duration + NeuralBuffer;
[~, ~, ~,~, ~, ~, ~, ISI] = process_SAT2(ss,NCells, SpikesArrivalTimes_Behav, SpikesArrivalTimes_Baseline, 'Operant','self',TimeOn, TimeOff, Who, ExpType,What);
FIG=figure();
H=histogram(ISI.Trill,0:5:1110, 'Normalization','probability');
ISI_Trill_df.time = H.BinEdges(1:end-1);
ISI_Trill_df.value = H.Values;
H=histogram(ISI.NonTrill,0:5:1110, 'Normalization','probability');
ISI_NonTrill_df.time = H.BinEdges(1:end-1);
ISI_NonTrill_df.value = H.Values;
close(FIG)

% Now go through all the recording during operant and highlight the
% probability of being in one or the other regime
Spike_arrival_times_ms = (Spike_arrival_times(Spike_arrival_times<OperantSession(2)*10^6) - OperantSession(1)*10^6).*10^-3; % alligned spike onset to start of the operant session and convert to ms
spike_regime(Spike_arrival_times_ms, 'Operant', NeuralWin, Overlap, Buffer, VocDataTime, ExpType, VocRank, What, ISI_Trill_df, ISI_NonTrill_df)

% Now go through all the recording during free session and highlight the
% probability of being in one or the other regime
Spike_arrival_times_ms = (Spike_arrival_times(Spike_arrival_times<FreeBehavSession(2)*10^6) - FreeBehavSession(1)*10^6).*10^-3; % alligned spike onset to start of the free session and convert to ms
spike_regime(Spike_arrival_times_ms, 'Free', NeuralWin, Overlap, Buffer, VocDataTime, ExpType, VocRank, What, ISI_Trill_df, ISI_NonTrill_df)

% Calculate the shape coefficent of the gamma distribution that best fits
% the distribution of ISI during operant


%% Explore the regimes of each cell
NeuralWin = 2000; %duration of the time window in ms
NeuralBuffer = 5000; %used to extract neural data by neuralData_compile_per_file when called by wrapper_allbehav_neuro_script
%       ms that should be used to calculate the time varying rate and
%       coefficient of variation
Overlap = 500; %ms
% Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
TR=1;% time resolution of Gaussian
nStd = 6; % before set as 4
Tau = (TR/2);
T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
Expwav = Expwav./sum(Expwav);
Buffer = 500; % if a vocalization is within Buffer ms of the edges of the NeuralWin window, then it is counted as during a vocalization
NCells = length(GoodCellIndices);
% Operant.Rate_CVISI_All = cell(NCells,3);
% Operant.Rate_CVISI_Voc = cell(NCells,3);
% Operant.Rate_CVISI_NonVoc = cell(NCells,3);
% Operant.Rate_CVISI_Trill = cell(NCells,3);
% Operant.Rate_CVISI_Ba = cell(NCells,3);
% Operant.SpikeSpectrum = cell(NCells,2);
% CellName_list = cell(NCells,1);
% Operant.Spectrum_Win_Voc = cell(NCells,2);
% Operant.Spectrum_Win_NonVoc = cell(NCells,2);
% Operant.Spectrum_All = cell(NCells,2);
% 
% FreeBehav = Operant;


for ss=211:NCells
    [Path2Cell,CellName] = fileparts(ListSSU{GoodCellIndices(ss)});
    fprintf(1,'File %d/%d: %s \n',ss,NCells,CellName)
    CellName_list{ss} = CellName;
    Path2Cell = [BasePath Path2Cell(strfind(Path2Cell, 'LMC'):end)];
    CellName_full = fullfile(Path2Cell, [CellName '.mat']);
    if exist(CellName_full, "file")
        load(CellName_full, 'Spike_arrival_times'); % load spike arrival times in useconds
    else
        warning("Issue with loading file %s cell %d/%d", CellName_full, ss, NCells)
        continue
    end
    Ind_ = strfind(CellName, '_');
    SUBJ = CellName(1:(Ind_(1)-1));
    Date = CellName(Ind_(1)+ (1:8));
    Tetrode = CellName(Ind_(2)+ 3);
    SSQ = CellName(Ind_(4) - 1);
    CellID = CellName((Ind_(4) + 1) : end);
    FileID = fullfile(OutputPath, [sprintf('%s_%s_SS%s_%s-%s',SUBJ, Date, SSQ, Tetrode, CellID) '.mat']);
    if exist(FileID,"file")
        load(FileID,'FreeBehavSession','OperantSession','PlayBackSession','QualitySSU','What','Who','ExpType','VocRank','Duration', 'DelayBefore','DelayAfter', 'SpikesArrivalTimes_Behav', 'SpikesArrivalTimes_Baseline') % these spike times are centered around vocalization onset
    else
        warning("Issue with loading file %s cell %d/%d", FileID, ss, NCells)
        continue
    end
    
    % Find vocalizations onset/offsets
    PathParts = regexp(Path2Cell, '/', 'split');
    Loggers_dir = ['/' fullfile(PathParts{1:end-2})];
    VocDataFiles= dir(fullfile(Loggers_dir, sprintf('%s_*_VocExtractData.mat', Date(3:end))));
    VocDataTime = cell(1,length(VocDataFiles));
    for ff=1:length(VocDataFiles)
        load(fullfile(VocDataFiles(ff).folder, VocDataFiles(ff).name), 'Voc_transc_time_refined'); % this variable is in ms
        % Make sure to re-order values in chronological order
        [~,SortInd] = sort(Voc_transc_time_refined(:,1));
        VocDataTime{ff} = Voc_transc_time_refined(SortInd,:)';
    end
    VocDataTime = [VocDataTime{:}]'- OperantSession(1)*10^3; % Operant Session value is in s and VocDataTime in ms
    
    % Fist look at spike regime for all the recording during operant
    % conditioning
    Spike_arrival_times_ms = (Spike_arrival_times(Spike_arrival_times<OperantSession(2)*10^6) - OperantSession(1)*10^6).*10^-3; % alligned spike onset to start of the operant session and convert to ms
    if length(Spike_arrival_times_ms)>1000
        % [Operant.Spectrum_Win_Voc(ss,:), Operant.Spectrum_Win_NonVoc(ss,:), Operant.Spectrum_All{ss,1}, Operant.Spectrum_All{ss,2}, Operant.Rate_CVISI_All(ss,:), Operant.Rate_CVISI_Voc(ss,:), Operant.Rate_CVISI_NonVoc(ss,:), Operant.Rate_CVISI_Trill(ss,:), Operant.Rate_CVISI_Ba(ss,:)] = process_SAT(ss,NCells, Spike_arrival_times_ms, 'operant', NeuralWin, Overlap, Buffer,VocDataTime, VocRank, What);
        TimeOn = -NeuralBuffer;
        TimeOff = Duration + NeuralBuffer;
        [Operant.Spectrum_Win_Voc(ss,:), Operant.Spectrum_Win_NonVoc(ss,:), Operant.Spectrum_Trill(ss,:),Operant.Spectrum_NonTrill(ss,:), Operant.Rate_CVISI_NonVoc(ss,:), Operant.Rate_CVISI_Voc(ss,:), Operant.Rate_CVISI_Trill(ss,:), Operant.ISI] = process_SAT2(ss,NCells, SpikesArrivalTimes_Behav, SpikesArrivalTimes_Baseline, 'Operant','self',TimeOn, TimeOff, Who, ExpType,What);
        
    else
        fprintf(1, 'Not enough spikes (min 1000) to run on operant data\n')
    end

    % Second look at spike regime for all the recording during Free session
    Spike_arrival_times_ms = (Spike_arrival_times(Spike_arrival_times<FreeBehavSession(2)*10^6) - FreeBehavSession(1)*10^6).*10^-3; % alligned spike onset to start of the free behav session and convert to ms
    if length(Spike_arrival_times_ms)>1000
        % [FreeBehav.Spectrum_Win_Voc(ss,:), FreeBehav.Spectrum_Win_NonVoc(ss,:), FreeBehav.Spectrum_All{ss,1}, FreeBehav.Spectrum_All{ss,2}, FreeBehav.Rate_CVISI_All(ss,:), FreeBehav.Rate_CVISI_Voc(ss,:), FreeBehav.Rate_CVISI_NonVoc(ss,:), FreeBehav.Rate_CVISI_Trill(ss,:), FreeBehav.Rate_CVISI_Ba(ss,:)] = process_SAT(ss, NCells, Spike_arrival_times_ms, 'free', NeuralWin, Overlap, Buffer,VocDataTime, VocRank, What);
        TimeOn = -NeuralBuffer;
        TimeOff = Duration + NeuralBuffer;
        [FreeBehav.Spectrum_Win_Voc(ss,:), FreeBehav.Spectrum_Win_NonVoc(ss,:), FreeBehav.Spectrum_Trill(ss,:), FreeBehav.Spectrum_NonTrill(ss,:), FreeBehav.Rate_CVISI_NonVoc(ss,:), FreeBehav.Rate_CVISI_Voc(ss,:), FreeBehav.Rate_CVISI_Trill(ss,:)] = process_SAT2(ss,NCells, SpikesArrivalTimes_Behav, SpikesArrivalTimes_Baseline, 'Free','self',TimeOn, TimeOff, Who, ExpType,What);
    else
        fprintf(1, 'Not enough spikes (min 1000) to run on free behavior data\n')
    end
        
    pause(1)
    
end
fprintf(' DONE \n')
save(fullfile(OutputPath, 'CellSpikingRegime2024.mat'),'ListSSU', 'GoodCellIndices', 'CellName_list', 'NeuralWin','Overlap', 'Operant','FreeBehav');


%% Run a umap on the CV_Rate plots
CapCV = 36;
RateCVISI_mat = nan(size(Rate_CVISI_All,1), numel(Rate_CVISI_All{1,1}));
for cc=1:NCells
    RateCVISI_mat(cc,:) = reshape(Rate_CVISI_All{cc,1}, 1, numel(Rate_CVISI_All{cc,1}))./sum(sum(Rate_CVISI_All{cc,1}));
end
% Finding the good cell indices (some cells have no data)
GoodCells = find(~isnan(RateCVISI_mat(:,1)));
figure(4)
clf
imagesc(RateCVISI_mat(GoodCells,:))
colorbar()
xlabel('CV on ISI and Rate profile')
ylabel('Cell#')

% Project the data in UMAP space
figure(5)
clf
[Reduction,UMAP,KMeans_ID]= run_umap(RateCVISI_mat(GoodCells,:), 'n_neighbors',9,'min_dist',0.051);

% Cluster the projected data with HDBSCAN
Cluster = HDBSCAN(Reduction);
Cluster.run_hdbscan(5,5,[],0.9);

figure(11)
clf
subplot(1,2,1)
Cluster.plot_tree()
subplot(1,2,2)
Cluster.plot_clusters()
suplabel('HDBSCAN clustering','t')


Yval = Rate_CVISI_All{1,3}(1:end-1) + diff(Rate_CVISI_All{1,3});
Xval = Rate_CVISI_All{1,2}(1:end-1) + diff(Rate_CVISI_All{1,2});

UkID = unique(KMeans_ID);
Cval = linspace(1,100,length(UkID));
Cval = Cval(randperm(length(UkID)));
Legend = cell(length(UkID),1);
Nr = floor(length(UkID)^0.5);
Nc = ceil(length(UkID)/Nr);
figure(6)
clf
figure(7)
clf
for cc=1:length(UkID)
   AvMap = mean(RateCVISI_mat(GoodCells(KMeans_ID==UkID(cc)),:));
   AvMap = reshape(AvMap, size(Rate_CVISI_All{1,1}))';
   figure(6)
   subplot(Nr,Nc,cc)
   imagesc(Xval, Yval, AvMap)
   colormap("default")
   axis xy
   ylabel('CV_ISI')
   set(gca, 'YLim',[0 2])
   XLim = get(gca, 'XLim');
   set(gca, 'XLim', [0 XLim(2)])
   set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
   xlabel('Spike Rate (Hz)')
   title(sprintf('Cluster %d (n=%d)',UkID(cc), sum(KMeans_ID==UkID(cc)) )) 
    
   figure(7)
   hold on
   scatter(Reduction(KMeans_ID==UkID(cc),1), Reduction(KMeans_ID==UkID(cc),2), 30, Cval(cc)*ones(length(Reduction(KMeans_ID==UkID(cc),1)),1), 'filled')
   Legend{cc} = (sprintf('Cluster %d (n=%d)',UkID(cc), sum(KMeans_ID==UkID(cc))));
end
figure(7)
xlabel('UMAP D1')
ylabel('UMAP D2')
% scatter(0, 0, 30, Cval(cc+3), 'filled')
% legend([Legend; ' '], 'Location','eastoutside')
legend(Legend, 'Location','eastoutside')
colormap(jet)
title('KMeans clustering')
hold off
figure(6)
suplabel('KMeans Average cluster map','t')


UHID = unique(Cluster.labels);
Cval = linspace(1,100,length(UHID));
Cval = Cval(randperm(length(UHID)));
Legend = cell(length(UHID),1);
Nr = floor(length(UHID)^0.5);
Nc = ceil(length(UHID)/Nr);
figure(8)
clf
figure(9)
clf

for cc=1:length(UHID)
   AvMap = mean(RateCVISI_mat(GoodCells(Cluster.labels==UHID(cc)),:));
   AvMap = reshape(AvMap, size(Rate_CVISI_All{1,1}))';
   figure(8)
   subplot(Nr,Nc,cc)
   imagesc(Xval, Yval, AvMap)
   colormap("default")
   axis xy
   ylabel('CV_ISI')
   set(gca, 'YLim',[0 2])
   XLim = get(gca, 'XLim');
   set(gca, 'XLim', [0 XLim(2)])
   set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
   xlabel('Spike Rate (Hz)')
   title(sprintf('Cluster %d (n=%d)',UHID(cc),sum(Cluster.labels==UHID(cc))))
   
   figure(9)
   hold on
   scatter(Reduction(Cluster.labels==UHID(cc),1), Reduction(Cluster.labels==UHID(cc),2), 30, Cval(cc)*ones(length(Reduction(Cluster.labels==UHID(cc),1)),1), 'filled')
   Legend{cc} = (sprintf('Cluster %d (n=%d)',UHID(cc),sum(Cluster.labels==UHID(cc))));
end
figure(9)
xlabel('UMAP D1')
ylabel('UMAP D2')
% scatter(0, 0, 30, Cval(cc+3), 'filled')
% legend([Legend; ' '], 'Location','eastoutside')
legend(Legend, 'Location','eastoutside')
colormap(jet)
title('HDBSCAN clustering')
hold off
   
figure(8)   
suplabel('HDBSCAN Average cluster map','t')

save(fullfile(OutputPath, 'CellSpikingRegime.mat'),'NCells', 'GoodCells','RateCVISI_mat','Reduction','UMAP','KMeans_ID','Cluster','-append');
%% OLD
figure()
subplot(1,3,1)
Yval = Rate_CVISI_All{122,3}(1:end-1) + diff(Rate_CVISI_All{122,3});
Xval = Rate_CVISI_All{122,2}(1:end-1) + diff(Rate_CVISI_All{122,2});
imagesc(Xval,Yval,Rate_CVISI_All{122,1}')
axis xy
ylabel('CV_ISI')
   XLim = get(gca, 'XLim');
   set(gca, 'XLim', [0 XLim(2)])
   set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
   xlabel('Spike Rate (Hz)')
   colorbar()
   
   subplot(1,3,2)
Yval = Rate_CVISI_All{122,3}(1:end-1) + diff(Rate_CVISI_All{122,3});
Xval = Rate_CVISI_All{122,2}(1:end-1) + diff(Rate_CVISI_All{122,2});
imagesc(Xval,Yval,Rate_CVISI_All{122,1}'./sum(sum(Rate_CVISI_All{122,1})))
axis xy
ylabel('CV_ISI')
   XLim = get(gca, 'XLim');
   set(gca, 'XLim', [0 XLim(2)])
   set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
   xlabel('Spike Rate (Hz)')
   colorbar()
   
   subplot(1,3,3)
imagesc(Xval, Yval,reshape(RateCVISI_mat(122,:), size(Rate_CVISI_All{122,1}))')
axis xy
   ylabel('CV_ISI')
   XLim = get(gca, 'XLim');
   set(gca, 'XLim', [0 XLim(2)])
   set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
   xlabel('Spike Rate (Hz)')
   colorbar()
%%   INTERNAL FUNCTIONS 
function [ListSSU] = gather_neural_datapath(BasePath)
fprintf(1,'*** Gathering paths to spike sorted units ***')
ListSSU = cell(10^3,1); % initialize the list to 1000
ExpFolders = dir(fullfile(BasePath,'LMC*'));
NSSU = 0; % counter for single units
for ee=1:length(ExpFolders)
    fprintf(1, '\n  -> Looking into  %s...\n ', fullfile(ExpFolders(ee).folder,ExpFolders(ee).name))
    DateFolders = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name, 'logger','20*'));
    for dd=1:length(DateFolders)
        fprintf(1, '   %s\n', DateFolders(dd).name);
        LoggerFolders = dir(fullfile(DateFolders(dd).folder, DateFolders(dd).name,'Logger*'));
        for ll = 1:length(LoggerFolders)
            SSFiles = dir(fullfile(LoggerFolders(ll).folder, LoggerFolders(ll).name, 'extracted_data', '*_TT*_SSS*.mat'));
            if ~isempty(SSFiles)
                for ssf=1:length(SSFiles)
                    NSSU = NSSU +1;
                    ListSSU{NSSU} = fullfile(SSFiles(ssf).folder, SSFiles(ssf).name);
                end
            end
            SSFiles = dir(fullfile(LoggerFolders(ll).folder, LoggerFolders(ll).name, 'extracted_data', '*_TT*_SSM*.mat'));
            if ~isempty(SSFiles)
                for ssf=1:length(SSFiles)
                    NSSU = NSSU +1;
                    ListSSU{NSSU} = fullfile(SSFiles(ssf).folder, SSFiles(ssf).name);
                end
            end
            SSFiles = dir(fullfile(LoggerFolders(ll).folder, LoggerFolders(ll).name, 'extracted_data', '*_TT*_SSU*.mat'));
            if ~isempty(SSFiles)
                for ssf=1:length(SSFiles)
                    NSSU = NSSU +1;
                    ListSSU{NSSU} = fullfile(SSFiles(ssf).folder, SSFiles(ssf).name);
                end
            end
        end
    end
end
ListSSU = ListSSU(1:NSSU);
fprintf(1, '\n Files from %d single units have been retrieved\n', NSSU);
end

function [] = spike_regime(Spike_arrival_times_ms, Session, NeuralWin, Overlap, Buffer, VocDataTime, ExpType, VocRank, What, ISI_Trill_df, ISI_NonTrill_df)
    FigTag=0; %Set to 1 to see figure of ISI fit for each windo    
    ISI_ms = diff(Spike_arrival_times_ms);
    TimeOn = -2*NeuralWin:Overlap:(max(Spike_arrival_times_ms)-NeuralWin);
    TimeOff = -1*NeuralWin:Overlap:max(Spike_arrival_times_ms);
    TR = NeuralWin/10; % time resolution to estimate the fano factor
    Rate = nan(length(TimeOn),1);
    Voc01 = nan(length(TimeOn),1);
    Regime_Trill_p = nan(length(TimeOn),1);
    Regime_NonTrill_p = nan(length(TimeOn),1);
    CoeffVar = nan(length(TimeOn),1);
    FanoFactor = nan(length(TimeOn),1);
    GammaA = nan(length(TimeOn),1);
    SSEGam = nan(length(TimeOn),1);
    SSEExp = nan(length(TimeOn),1);
    LLGam = nan(length(TimeOn),2);
    LLExp = nan(length(TimeOn),2);
    LLRT.pValue = nan(length(TimeOn),1);
    LLRT.stat = nan(length(TimeOn),1);
    if strcmp(Session, 'Operant') && sum(contains(VocRank(contains(ExpType, 'O')), 'first')) == size(VocDataTime,1) % we can retrieve the type of call!
        Trill01 = zeros(length(TimeOn),1); % Trill =1; Non-Trill or non Voc = 0 -> Non-Trill = Voc01 .* ~Trill01
        What_local = What(logical(contains(ExpType, 'O') .* contains(VocRank, 'first')));
    end
    for bb=1:length(TimeOn)
        Spikes01 = (Spike_arrival_times_ms>=TimeOn(bb)) .* (Spike_arrival_times_ms<TimeOff(bb));
        Rate(bb) = sum(Spikes01)/(NeuralWin*10^-3);
        ISI_local = ISI_ms(logical(Spikes01(1:(end-1))));
        Spike_arrival_times_ms_local = Spike_arrival_times_ms(logical(Spikes01))-TimeOn(bb);
        % get the probability of obtaining these ISI under Regime Trill and
        % Regime NonTrill
        if isempty(ISI_local)
            Regime_Trill_p(bb) = 0;
            Regime_NonTrill_p(bb) = 0;
        else
            RT_p_local = nan(length(ISI_local),1);
            RNT_p_local = nan(length(ISI_local),1);
            for ii=1:length(ISI_local)
                IndexTrill = find(ISI_Trill_df.time>=ISI_local(ii),1)-1;
                if isempty(IndexTrill) % the value is beyond the distribution definition
                    RT_p_local(ii) = ISI_Trill_df.value(end)/2;
                elseif IndexTrill==0
                    RT_p_local(ii) = ISI_Trill_df.value(1);
                else
                    RT_p_local(ii) = ISI_Trill_df.value(IndexTrill);
                end
                IndexNonTrill = find(ISI_NonTrill_df.time>=ISI_local(ii),1)-1;
                if isempty(IndexNonTrill) % the value is beyond the distribution definition
                    RNT_p_local(ii) = ISI_NonTrill_df.value(end)/2;
                elseif IndexNonTrill==0
                    RNT_p_local(ii) = ISI_NonTrill_df.value(1);
                else
                    RNT_p_local(ii) = ISI_NonTrill_df.value(IndexNonTrill);
                end
            end
            Regime_Trill_p(bb) = sum(RT_p_local)./length(RT_p_local);
            Regime_NonTrill_p(bb) = sum(RNT_p_local)./length(RNT_p_local);
        end
        
        % Calculate and plot estimates of the noise distribution (Regime
        % Poisson or Gamma)
        if ~isempty(ISI_local) && (length(ISI_local)>=10)
            % This is the coefficient of variation of ISI
            CoeffVar(bb) = std(ISI_local)/mean(ISI_local);
            % This is the Fano Factor of the spike count per bin TR (=1 for a poisson process)
            TimeBins = 0:TR:NeuralWin;
            SpikeCountPerBin = nan(length(TimeBins)-1,1);
            for tt=1:(length(TimeBins)-1)
                SpikeCountPerBin(tt) = sum((Spike_arrival_times_ms_local>TimeBins(tt)) .* (Spike_arrival_times_ms_local<TimeBins(tt+1)));
            end
            FanoFactor(bb) = (std(SpikeCountPerBin).^2) ./ mean(SpikeCountPerBin);
            % plot the histogram of ISI
            if FigTag
                figure(20);subplot(2,1,1);cla
                H=histogram(ISI_local, 'BinWidth',1); xlabel('ISI (ms)');
            else
                [H.BinCounts,H.BinEdges] =histcounts(ISI_local, 'BinWidth',1);
            end
            % fit the distribution of ISI with a gamma distribution
            Gam = fitdist(ISI_local, 'Gamma');
            GammaA(bb) = Gam.a; % shape of the gamma distribution fitted
            % fit the distribution of ISI with an exponential (Poisson process)
            Poi = fitdist(ISI_local, 'Exponential');
            % Calculate the posterior probabilities of observing this sequence of ISI given the Poisson or the Gamma model and add the2 fit to the figure
            XISI = H.BinEdges(1):(H.BinEdges(end)-1);
            Xexp = pdf('Exponential',XISI,Poi.mu);
            Xgam = gampdf(XISI, Gam.a, Gam.b);
            pISI_local_exp = pdf('Exponential',ISI_local,Poi.mu);
            pISI_local_gam = gampdf(ISI_local, Gam.a, Gam.b);
            % Calculate SSE and log likelihood
            PdfH = H.BinCounts./sum(H.BinCounts);
            SSEGam(bb) = sum((PdfH-Xgam).^2);
            SSEExp(bb) = sum((PdfH-Xexp).^2);
            LLGam(bb,1) = sum(log(pISI_local_gam));
            LLGam(bb,2) = length(pISI_local_gam);
            LLExp(bb,1) = sum(log(pISI_local_exp));
            LLExp(bb,2) = length(pISI_local_exp); % keep track of the number of element per comparison to plot normalized LL
            LLRT.stat(bb) = -2*(LLExp(bb,1)-LLGam(bb,1));
            LLRT.pValue(bb) = chi2cdf(LLRT.stat(bb), 1, 'upper');
            
            % plot everything
            if FigTag
                figure(20);subplot(2,1,1); cla
                yyaxis right; plot(XISI, Xexp, '-', 'LineWidth',2, 'Color', [0.929 0.694 0.125]);
                hold on; yyaxis right; plot(XISI, Xgam, '-', 'LineWidth',2);
                hold on; yyaxis left; H=histogram(ISI_local, 'BinWidth',1); xlabel('ISI (ms)');hold off
                legend({'ISI','Exp fit', 'Gamma fit'})
                title(sprintf('ISI fit Gamma parameter = %.2f SSEGam = %.2e SSEExp = %.2e FanoFactor = %.2f', Gam.a, SSEGam(bb), SSEExp(bb), FanoFactor(bb)))
            end
        end
        % identify if a vocalization was produced during that time window
        if isempty(VocDataTime)
            Voc01(bb) = 0;
        else
            Voc01(bb) = any(((VocDataTime(:,1)-Buffer)<=TimeOn(bb)) .* ((VocDataTime(:,2)+Buffer)>=TimeOn(bb))) || any(((VocDataTime(:,1)-Buffer)<=TimeOff(bb)) .* ((VocDataTime(:,2)+ Buffer)>=TimeOff(bb))) || any(((VocDataTime(:,1)-Buffer)>=TimeOn(bb)) .* ((VocDataTime(:,2)+Buffer)<=TimeOff(bb)));
            if strcmp(Session, 'Operant') && Voc01(bb) && (sum(contains(VocRank(contains(ExpType, 'O')), 'first')) == size(VocDataTime,1)) % we can retrieve the type of call!
                VocInd = [find(((VocDataTime(:,1)-Buffer)<=TimeOn(bb)) .* ((VocDataTime(:,2)+Buffer)>=TimeOn(bb))) find(((VocDataTime(:,1)-Buffer)<=TimeOff(bb)) .* ((VocDataTime(:,2)+ Buffer)>=TimeOff(bb))) find(((VocDataTime(:,1)-Buffer)>=TimeOn(bb)) .* ((VocDataTime(:,2)+Buffer)<=TimeOff(bb)))];
                VocInd = unique(VocInd);
                if length(VocInd)>2
                    keyboard
                end
                Trill01(bb) = contains(What_local{VocInd}, 'Tr');
            end
        end
    end
    
    % remove points due to lost cell
    % Define deadtime for the cell as 1 min of rate =0
    DeadPoints = strfind(((Rate(1:end-1)==0).*(diff(Rate)==0))', ones(1,round((60*10^3)/NeuralWin)));
    if isempty(DeadPoints)
        FirstDead = length(Rate);
    else
        FirstDead = DeadPoints(1);
    end
    Rate = Rate(1:(FirstDead-2));
    Voc01 = Voc01(1:(FirstDead-2));
    if strcmp(Session, 'Operant') && sum(contains(VocRank(contains(ExpType, 'O')), 'first')) == size(VocDataTime,1)
        Trill01 = Trill01(1:(FirstDead-2));
    end
    TimeOn = TimeOn(1:FirstDead-2);
    Regime_Trill_p = Regime_Trill_p(1:FirstDead-2);
    Regime_NonTrill_p = Regime_NonTrill_p(1:FirstDead-2);
    Index0 = Regime_Trill_p==0;
    Regime_Trill_p(Index0) = min(Regime_Trill_p(~Index0));
    Index0 = Regime_NonTrill_p==0;
    Regime_NonTrill_p(Index0) = min(Regime_NonTrill_p(~Index0));
    CoeffVar = CoeffVar(1:FirstDead-2);
    FanoFactor = FanoFactor(1:FirstDead-2);
    GammaA = GammaA(1:FirstDead-2);
    SSEGam = SSEGam(1:FirstDead-2);
    SSEExp = SSEExp(1:FirstDead-2);
    LLGam = LLGam(1:FirstDead-2,:);
    LLExp = LLExp(1:FirstDead-2,:);
    LLRT.pValue = LLRT.pValue(1:FirstDead-2);
    LLRT.stat = LLRT.stat(1:FirstDead-2);

    figure(2)
    
    clf
    subplot(2,1,1)
    y = Rate';
    x=(TimeOn+NeuralWin/2).*10^-3;
    z = zeros(size(x));
    lineColor = log2(Regime_Trill_p ./ Regime_NonTrill_p);  % This is the color, it varies with x in this case.
    % Plot the line with width 8 so we can see the colors well.
    surface([x;x], [y;y], [z;z], [lineColor';lineColor'],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 2);
    xlabel('Time (s)')
    ylabel('Spike rate (Hz)')
    title(sprintf('%s session Bin size %ds step %ds', Session, NeuralWin/10^3, Overlap/10^3))
    colormap('cool')
    c=colorbar();
    c.Label.String = 'log ratio of probability: log2(Trill/NonTrill)';
    CL = max(abs(lineColor));
    clim([-CL CL])
    xlim([0 ceil(max((TimeOn+NeuralWin/2).*10^-3))])
    hold on
    if exist('Trill01','var')
        plot((TimeOn(logical(Trill01))+NeuralWin/2).*10^-3,ones(sum(Trill01),1), '*r')
        hold on
        plot((TimeOn(logical(Voc01 .* ~Trill01))+NeuralWin/2).*10^-3,ones(sum(Voc01 .* ~Trill01),1), '*b')
        hold off
        legend('Rate', 'Trill', 'NonTrill')
    else
        plot((TimeOn(logical(Voc01))+NeuralWin/2).*10^-3,ones(sum(Voc01),1), '*k')
        hold off
        legend('Rate', 'Voc')
    end
    
    subplot(2,1,2)
    hold on
    plot((TimeOn+NeuralWin/2).*10^-3, Regime_NonTrill_p, '-b','LineWidth',2 )
    plot((TimeOn+NeuralWin/2).*10^-3, Regime_Trill_p, '-r','LineWidth',2 )
    hold off
    xlabel('Time (s)')
    ylabel('Probability')
    legend('pRegime NonTrill', 'pRegime Trill')
    xlim([0 ceil(max((TimeOn+NeuralWin/2).*10^-3))])
    hold off
    
    figure(5)
    clf
    subplot(4,1,1)
    y = Rate';
    x=(TimeOn+NeuralWin/2).*10^-3;
    z = zeros(size(x));
    lineColor = CoeffVar;  % This is the color, it varies with x in this case.
    % Plot the line with width 8 so we can see the colors well.
    surface([x;x], [y;y], [z;z], [lineColor';lineColor'],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 2);
    xlabel('Time (s)')
    ylabel('Spike rate (Hz)')
    title(sprintf('%s session Bin size %ds step %ds', Session, NeuralWin/10^3, Overlap/10^3))
    colormap('cool')
    c=colorbar();
    c.Label.String = 'Coefficient of variation of ISI (std/mean)';
    xlim([0 ceil(max((TimeOn+NeuralWin/2).*10^-3))])
    hold on
    if exist('Trill01','var')
        plot((TimeOn(logical(Trill01))+NeuralWin/2).*10^-3,ones(sum(Trill01),1), '*r')
        hold on
        plot((TimeOn(logical(Voc01 .* ~Trill01))+NeuralWin/2).*10^-3,ones(sum(Voc01 .* ~Trill01),1), '*b')
        hold off
        legend('Rate', 'Trill', 'NonTrill')
    else
        plot((TimeOn(logical(Voc01))+NeuralWin/2).*10^-3,ones(sum(Voc01),1), '*k')
        hold off
        legend('Rate', 'Voc')
    end

    subplot(4,1,2)
    y = Rate';
    x=(TimeOn+NeuralWin/2).*10^-3;
    z = zeros(size(x));
    lineColor = log10(FanoFactor);  % This is the color, it varies with x in this case.
    % Plot the line with width 8 so we can see the colors well.
    surface([x;x], [y;y], [z;z], [lineColor';lineColor'],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 2);
    xlabel('Time (s)')
    ylabel('Spike rate (Hz)')
    title(sprintf('%s session Bin size %ds step %ds', Session, NeuralWin/10^3, Overlap/10^3))
    colormap('cool')
    c=colorbar();
    c.Label.String = 'log10 Fano Factor of spike count (var/mean)';
    CL = max(abs(lineColor));
    clim([-CL CL])
    xlim([0 ceil(max((TimeOn+NeuralWin/2).*10^-3))])
    hold on
    if exist('Trill01','var')
        plot((TimeOn(logical(Trill01))+NeuralWin/2).*10^-3,ones(sum(Trill01),1), '*r')
        hold on
        plot((TimeOn(logical(Voc01 .* ~Trill01))+NeuralWin/2).*10^-3,ones(sum(Voc01 .* ~Trill01),1), '*b')
        hold off
        legend('Rate', 'Trill', 'NonTrill')
    else
        plot((TimeOn(logical(Voc01))+NeuralWin/2).*10^-3,ones(sum(Voc01),1), '*k')
        hold off
        legend('Rate', 'Voc')
    end

    subplot(4,1,3)
    y = Rate';
    x=(TimeOn+NeuralWin/2).*10^-3;
    z = zeros(size(x));
    lineColor = GammaA;  % This is the color, it varies with x in this case.
    % Plot the line with width 8 so we can see the colors well.
    surface([x;x], [y;y], [z;z], [lineColor';lineColor'],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 2);
    xlabel('Time (s)')
    ylabel('Spike rate (Hz)')
    title(sprintf('%s session Bin size %ds step %ds', Session, NeuralWin/10^3, Overlap/10^3))
    colormap('cool')
    c=colorbar();
    c.Label.String = 'Shape of the gamma fit to ISI';
    % CL = max(abs(lineColor));
    clim([0 10])
    xlim([0 ceil(max((TimeOn+NeuralWin/2).*10^-3))])
    hold on
    if exist('Trill01','var')
        plot((TimeOn(logical(Trill01))+NeuralWin/2).*10^-3,ones(sum(Trill01),1), '*r')
        hold on
        plot((TimeOn(logical(Voc01 .* ~Trill01))+NeuralWin/2).*10^-3,ones(sum(Voc01 .* ~Trill01),1), '*b')
        hold off
        legend('Rate', 'Trill', 'NonTrill')
    else
        plot((TimeOn(logical(Voc01))+NeuralWin/2).*10^-3,ones(sum(Voc01),1), '*k')
        hold off
        legend('Rate', 'Voc')
    end

    subplot(4,1,4)
    y = Rate';
    x=(TimeOn+NeuralWin/2).*10^-3;
    z = zeros(size(x));
    SSEGam(isinf(SSEGam))=2*max(SSEGam(~isinf(SSEGam)));
    lineColor = log2(SSEExp ./ SSEGam);  % This is the color, it varies with x in this case.
    % Plot the line with width 8 so we can see the colors well.
    surface([x;x], [y;y], [z;z], [lineColor';lineColor'],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 2);
    xlabel('Time (s)')
    ylabel('Spike rate (Hz)')
    title(sprintf('%s session Bin size %ds step %ds', Session, NeuralWin/10^3, Overlap/10^3))
    colormap('cool')
    c=colorbar();
    c.Label.String = 'SSE Exp / Gamma';
    CL = max(abs(lineColor));
    clim([-CL CL])
    xlim([0 ceil(max((TimeOn+NeuralWin/2).*10^-3))])
    hold on
    if exist('Trill01','var')
        plot((TimeOn(logical(Trill01))+NeuralWin/2).*10^-3,ones(sum(Trill01),1), '*r')
        hold on
        plot((TimeOn(logical(Voc01 .* ~Trill01))+NeuralWin/2).*10^-3,ones(sum(Voc01 .* ~Trill01),1), '*b')
        hold off
        legend('Rate', 'Trill', 'NonTrill')
    else
        plot((TimeOn(logical(Voc01))+NeuralWin/2).*10^-3,ones(sum(Voc01),1), '*k')
        hold off
        legend('Rate', 'Voc')
    end
end

function [Spectrum_Win_Voc_out, Spectrum_Win_NonVoc_out, Spectrum, Freqs, Rate_CVISI_All, Rate_CVISI_Voc, Rate_CVISI_NonVoc, Rate_CVISI_Trill, Rate_CVISI_Ba] = process_SAT(ss, NCells, Spike_arrival_times_ms, SessionName, NeuralWin, Overlap, Buffer, VocDataTime, VocRank, What)
    ISI_ms = diff(Spike_arrival_times_ms);
    TimeOn = -2*NeuralWin:Overlap:(max(Spike_arrival_times_ms)-NeuralWin);
    TimeOff = -1*NeuralWin:Overlap:max(Spike_arrival_times_ms);
    Rate = nan(length(TimeOn),1);
    STD_ISI = nan(length(TimeOn),1);
    Mean_ISI = nan(length(TimeOn),1);
    Spectrum_Win = cell(length(TimeOn),2);
    Voc01 = nan(length(TimeOn),1);
    if sum(contains(VocRank, 'first')) == size(VocDataTime,1) % we can retrieve the type of call!
        Trill01 = zeros(length(TimeOn),1); % Trill =1; Non-Trill or non Voc = 0 -> Non-Trill = Voc01 .* ~Trill01
    end
    for bb=1:length(TimeOn)
        Spikes01 = (Spike_arrival_times_ms>=TimeOn(bb)) .* (Spike_arrival_times_ms<TimeOff(bb));
        [Spectrum_Win{bb,1}, Spectrum_Win{bb,2}] = spike_pattern_spectrum(Spike_arrival_times_ms(logical(Spikes01)), TimeOn(bb), TimeOff(bb));
        Rate(bb) = sum(Spikes01)/(NeuralWin*10^-3);
        STD_ISI(bb) = std(ISI_ms(logical(Spikes01(1:(end-1)))));
        Mean_ISI(bb) = mean(ISI_ms(logical(Spikes01(1:(end-1)))));
        if isempty(VocDataTime)
            Voc01(bb) = 0;
        else
            Voc01(bb) = any(((VocDataTime(:,1)-Buffer)<=TimeOn(bb)) .* ((VocDataTime(:,2)+Buffer)>=TimeOn(bb))) || any(((VocDataTime(:,1)-Buffer)<=TimeOff(bb)) .* ((VocDataTime(:,2)+ Buffer)>=TimeOff(bb))) || any(((VocDataTime(:,1)-Buffer)>=TimeOn(bb)) .* ((VocDataTime(:,2)+Buffer)<=TimeOff(bb)));
            if Voc01(bb) && (sum(contains(VocRank, 'first')) == size(VocDataTime,1)) % we can retrieve the type of call!
                VocInd = [find(((VocDataTime(:,1)-Buffer)<=TimeOn(bb)) .* ((VocDataTime(:,2)+Buffer)>=TimeOn(bb))) find(((VocDataTime(:,1)-Buffer)<=TimeOff(bb)) .* ((VocDataTime(:,2)+ Buffer)>=TimeOff(bb))) find(((VocDataTime(:,1)-Buffer)>=TimeOn(bb)) .* ((VocDataTime(:,2)+Buffer)<=TimeOff(bb)))];
                VocInd = unique(VocInd);
                if length(VocInd)>2
                    keyboard
                end
                Trill01(bb) = contains(What{VocInd}, 'Tr');
            end
        end
    end
    CV_ISI = STD_ISI ./ Mean_ISI; % Maybe this should be var/Mean to find which distributions are Poisson?
    CV_ISI(Rate==0) = 0; % Correct NaN values of CV when Rate =0;
    
    % remove points due to lost cell
    % Define deadtime for the cell as 1 min of rate =0
    DeadPoints = strfind(((Rate(1:end-1)==0).*(diff(Rate)==0))', ones(1,round((60*10^3)/NeuralWin)));
    if isempty(DeadPoints)
        FirstDead = length(Rate);
    else
        FirstDead = DeadPoints(1);
    end
    Rate = Rate(1:(FirstDead-2));
    STD_ISI = STD_ISI(1:(FirstDead-2));
    CV_ISI = CV_ISI(1:(FirstDead-2));
    Mean_ISI = Mean_ISI(1:(FirstDead-2));
    Spectrum_Win = Spectrum_Win(1:(FirstDead-2),:);
    Voc01 = Voc01(1:(FirstDead-2));
    if sum(contains(VocRank, 'first')) == size(VocDataTime,1)
        Trill01 = Trill01(1:(FirstDead-2));
    end
    Spectrum_Win_Voc = Spectrum_Win(logical(Voc01),:);
    Spectrum_Win_NonVoc = Spectrum_Win(~Voc01,:);

    % Run an HDBSCAN on the Rate and CV data to find the different
    % "states" of the cell
    LogRate = Rate;
    LogRate(LogRate == 0) = 0.1; % Estimate the minimum rate at 0.1 spike per second
    LogRate = log10(LogRate);
    % Clusterer = HDBSCAN([LogRate CV_ISI]);
    % Clusterer.run_hdbscan(5,5,[],0.9);
    
    % Construct a time varying vector of the neural activity and calculate
    % its spectrum for the whole session
    [Spectrum, Freqs] = spike_pattern_spectrum(Spike_arrival_times_ms);


    Rate_CVISI_All = cell(1,3);
    Rate_CVISI_Voc = cell(1,3);
    Rate_CVISI_NonVoc = cell(1,3);
    Rate_CVISI_Trill = cell(1,3);
    Rate_CVISI_Ba = cell(1,3);

    if sum(contains(VocRank, 'first')) == size(VocDataTime,1)
        if strcmp(SessionName, 'operant')
            figure(1)
        else
            figure(4)
        end
        clf
        H=histogram2(LogRate, CV_ISI, 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5], 'NumBins',[24 36]);
        XBinEdges = H.XBinEdges;
        YBinEdges = H.YBinEdges;
        Rate_CVISI_All{1} = H.Values;
        Rate_CVISI_All{2} = H.XBinEdges;
        Rate_CVISI_All{3} = H.YBinEdges;
%         set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
%         xlabel('Spike Rate (Hz)')
%         ylabel('CV of ISI')
%         title('all time points')

        subplot(2,2,1)
        H=histogram2(LogRate(logical(Voc01)), CV_ISI(logical(Voc01)), 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5], 'XBinEdges',XBinEdges, 'YBinEdges',YBinEdges);
        Rate_CVISI_Voc{1} = H.Values;
        Rate_CVISI_Voc{2} = H.XBinEdges;
        Rate_CVISI_Voc{3} = H.YBinEdges;
        set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title(sprintf('Time points around vocalizations (within %d ms)', Buffer+Overlap))

        subplot(2,2,2)
        H=histogram2(LogRate(~Voc01), CV_ISI(~Voc01), 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5], 'XBinEdges',XBinEdges, 'YBinEdges',YBinEdges);
        Rate_CVISI_NonVoc{1} = H.Values;
        Rate_CVISI_NonVoc{2} = H.XBinEdges;
        Rate_CVISI_NonVoc{3} = H.YBinEdges;
        set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title('Time points outside of vocalizations')
        
        subplot(2,2,3)
        H=histogram2(LogRate(logical(Voc01.*Trill01)), CV_ISI(logical(Voc01.*Trill01)), 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5], 'XBinEdges',XBinEdges, 'YBinEdges',YBinEdges);
        Rate_CVISI_Trill{1} = H.Values;
        Rate_CVISI_Trill{2} = H.XBinEdges;
        Rate_CVISI_Trill{3} = H.YBinEdges;
        set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title(sprintf('Time points around Trills (within %d ms)', Buffer+Overlap))
        
        subplot(2,2,4)
        H=histogram2(LogRate(logical(Voc01.*~Trill01)), CV_ISI(logical(Voc01.*~Trill01)), 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5], 'XBinEdges',XBinEdges, 'YBinEdges',YBinEdges);
        Rate_CVISI_Ba{1} = H.Values;
        Rate_CVISI_Ba{2} = H.XBinEdges;
        Rate_CVISI_Ba{3} = H.YBinEdges;
        set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title(sprintf('Time points around Barks (within %d ms)', Buffer+Overlap))
    else
        if strcmp(SessionName, 'operant')
            figure(1)
        else
            figure(4)
        end
        clf
        subplot(1,3,1)
        H=histogram2(LogRate, CV_ISI, 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5],'NumBins',[24 36]);
        Rate_CVISI_All{1} = H.Values;
        Rate_CVISI_All{2} = H.XBinEdges;
        Rate_CVISI_All{3} = H.YBinEdges;
        set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title('all time points')

        subplot(1,3,2)
        H=histogram2(LogRate(logical(Voc01)), CV_ISI(logical(Voc01)), 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5], 'XBinEdges',H.XBinEdges, 'YBinEdges',H.YBinEdges);
        Rate_CVISI_Voc{1} = H.Values;
        Rate_CVISI_Voc{2} = H.XBinEdges;
        Rate_CVISI_Voc{3} = H.YBinEdges;
        set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title(sprintf('Time points around vocalizations (within %d ms)', Buffer+Overlap))

        subplot(1,3,3)
        H=histogram2(LogRate(~Voc01), CV_ISI(~Voc01), 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5], 'XBinEdges',H.XBinEdges, 'YBinEdges',H.YBinEdges);
        Rate_CVISI_NonVoc{1} = H.Values;
        Rate_CVISI_NonVoc{2} = H.XBinEdges;
        Rate_CVISI_NonVoc{3} = H.YBinEdges;
        set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title('Time points outside of vocalizations')
        suplabel(sprintf('%d/%d %s', ss, NCells,SessionName));
    end
    
    
    if strcmp(SessionName, 'operant')
        figure(2)
    else
        figure(5)
    end
    clf
    subplot(2,1,1)
    plot((TimeOn(1:FirstDead-2)+NeuralWin/2).*10^-3,Rate', '-k', 'LineWidth',2)
    xlabel('Time (s)')
    zlabel('Hz')
    title('Spike Rate (Hz)')
    % colorbar()
    subplot(2,1,2)
    plot((TimeOn(1:FirstDead-2)+NeuralWin/2).*10^-3,CV_ISI', '-k', 'LineWidth',2)
    xlabel('Time (s)')
    title('CV ISI')
    % colorbar()
    suplabel(sprintf('%d/%d %s', ss, NCells,SessionName));
    
    if strcmp(SessionName, 'operant')
        figure(3)
    else
        figure(6)
    end
    Spectrum_Win_Voc_out = cell(2,1);
    Spectrum_Win_NonVoc_out = cell(2,1);
    if ~isempty(Spectrum_Win_Voc)
        Spectrum_Win_Voc_out{2} = Spectrum_Win_Voc{1,2};
        Spectrum_Win_Voc_out{1} = mean(cell2mat(Spectrum_Win_Voc(:,1)));
    end
    if ~isempty(Spectrum_Win_NonVoc)
        Spectrum_Win_NonVoc_out{2} = Spectrum_Win_NonVoc{1,2};
        Spectrum_Win_NonVoc_out{1} = mean(cell2mat(Spectrum_Win_NonVoc(:,1)));
    end

    clf
    plot(Freqs(Freqs<=30),10*log10(Spectrum(Freqs<=30)), '-k', 'LineWidth',2)
    hold on
    if ~isempty(Spectrum_Win_Voc)
        plot(Spectrum_Win_Voc{1,2}, 10*log10(mean(cell2mat(Spectrum_Win_Voc(:,1)))),  '-','Color',[0.8500 0.3250 0.0980], 'LineWidth',2)
        hold on
    end
    if ~isempty(Spectrum_Win_NonVoc)
        plot(Spectrum_Win_NonVoc{1,2}, 10*log10(mean(cell2mat(Spectrum_Win_NonVoc(:,1)))), '-', 'Color', [0 0.4470 0.7410], 'LineWidth',2)
    end
    xlim([0 30]);
    hold off
    ylabel('Power of the Spike pattern (DB/Hz)')
    xlabel('Frequency (Hz)')
    title(sprintf('%d/%d Spike Power Spectrum %s', ss, NCells, SessionName))
    if ~isempty(Spectrum_Win_Voc) & ~isempty(Spectrum_Win_NonVoc)
        legend('all','Voc 2s window', 'NonVoc 2s window')
    elseif isempty(Spectrum_Win_Voc) & ~isempty(Spectrum_Win_NonVoc)
        legend('all', 'NonVoc 2s window')
    elseif ~isempty(Spectrum_Win_Voc) & isempty(Spectrum_Win_NonVoc)
        legend('all','Voc 2s window')
    elseif isempty(Spectrum_Win_Voc) & isempty(Spectrum_Win_NonVoc)
        legend('all')
    end
end

function [Spectrum_Win_Voc_out, Spectrum_Win_baseline_out, Spectrum_Win_Trill_out,Spectrum_Win_NonTrill_out, Rate_CVISI_baseline, Rate_CVISI_Voc, Rate_CVISI_Trill, ISI_ms] = process_SAT2(ss,NCells, SpikesArrivalTimes_Behav, SpikesArrivalTimes_Baseline, SessionName,ID, TimeOn, TimeOff, Who, ExpType,What)
    if strcmp(ID,'self')
        EventIndices = find(contains(ExpType, SessionName(1)) .* contains(Who, 'self').* contains(What, 'Voc'));
    else
        EventIndices = find(contains(ExpType, SessionName(1)) .* ~contains(Who, 'self').* contains(What, 'Voc'));
    end
    Nevents = length(EventIndices);
    Rate_Voc = nan(Nevents,1);
    STD_ISI_Voc = nan(Nevents,1);
    Mean_ISI_Voc = nan(Nevents,1);
    ISI_ms_Behav = cell(1,Nevents);
    ISI_ms_baseline = cell(1,Nevents);
    Rate_baseline = nan(Nevents,1);
    STD_ISI_baseline = nan(Nevents,1);
    Mean_ISI_baseline = nan(Nevents,1);
    Spectrum_Win_Voc = cell(Nevents,2);
    Spectrum_Win_Baseline = cell(Nevents,2);
    Trill01 = nan(Nevents,1);
    for ee=1:Nevents
        eei = EventIndices(ee);
        [Spectrum_Win_Voc{ee,1}, Spectrum_Win_Voc{ee,2}] = spike_pattern_spectrum(SpikesArrivalTimes_Behav{eei}, TimeOn, TimeOff(eei));
        [Spectrum_Win_Baseline{ee,1}, Spectrum_Win_Baseline{ee,2}] = spike_pattern_spectrum(SpikesArrivalTimes_Baseline{eei}, -6000, -5000);% fixed one second
        Rate_Voc(ee) = length(SpikesArrivalTimes_Behav{eei})/((TimeOff(eei)-TimeOn)*10^-3);
        Rate_baseline(ee) = length(SpikesArrivalTimes_Baseline{eei}); % this is the number of spike in a one sec window
        ISI_ms_Behav{ee} = diff(SpikesArrivalTimes_Behav{eei});
        STD_ISI_Voc(ee) = std(ISI_ms_Behav{ee});
        Mean_ISI_Voc(ee) = mean(ISI_ms_Behav{ee});
        ISI_ms_baseline{ee} = diff(SpikesArrivalTimes_Baseline{eei}');
        STD_ISI_baseline(ee) = std(ISI_ms_baseline{ee});
        Mean_ISI_baseline(ee) = mean(ISI_ms_baseline{ee});
        Trill01(ee) = contains(What{eei}, 'Tr');
    end
    CV_ISI_Voc = STD_ISI_Voc ./ Mean_ISI_Voc; % Maybe this should be var/Mean to find which distributions are Poisson?
    CV_ISI_Voc(Rate_Voc==0) = 0; % Correct NaN values of CV when Rate =0;
    CV_ISI_baseline = STD_ISI_baseline ./ Mean_ISI_baseline; % Maybe this should be var/Mean to find which distributions are Poisson?
    CV_ISI_baseline(Rate_baseline==0) = 0; % Correct NaN values of CV when Rate =0;
    
    % remove points due to lost cell
    % Define deadtime for the cell as 3 sequences with 3 rates in a row
    DeadPoints = strfind(((Rate_Voc(1:end-1)==0).*(diff(Rate_Voc)==0))', ones(1,3));
    if isempty(DeadPoints)
        FirstDead = length(Rate_Voc);
    else
        FirstDead = DeadPoints(1);
    end
    Rate_Voc = Rate_Voc(1:(FirstDead));
    STD_ISI_Voc = STD_ISI_Voc(1:(FirstDead));
    CV_ISI_Voc = CV_ISI_Voc(1:(FirstDead));
    Mean_ISI_Voc = Mean_ISI_Voc(1:(FirstDead));
    ISI_ms_Behav = ISI_ms_Behav(1:FirstDead);
    Spectrum_Win_Voc = Spectrum_Win_Voc(1:(FirstDead),:);
    Rate_baseline = Rate_baseline(1:(FirstDead));
    STD_ISI_baseline = STD_ISI_baseline(1:(FirstDead));
    CV_ISI_baseline = CV_ISI_baseline(1:(FirstDead));
    Mean_ISI_baseline = Mean_ISI_baseline(1:(FirstDead));
    ISI_ms_baseline = ISI_ms_baseline(1:FirstDead);
    Spectrum_Win_Baseline = Spectrum_Win_Baseline(1:(FirstDead),:);
    Trill01 = Trill01(1:(FirstDead));
    

    % Run an HDBSCAN on the Rate and CV data to find the different
    % "states" of the cell
    LogRate_Voc = Rate_Voc;
    LogRate_Voc(LogRate_Voc == 0) = 0.1; % Estimate the minimum rate at 0.1 spike per second
    LogRate_Voc = log10(LogRate_Voc);
    LogRate_baseline = Rate_baseline;
    LogRate_baseline(LogRate_baseline == 0) = 0.1; % Estimate the minimum rate at 0.1 spike per second
    LogRate_baseline = log10(LogRate_baseline);
    % Clusterer = HDBSCAN([LogRate CV_ISI]);
    % Clusterer.run_hdbscan(5,5,[],0.9);


    Rate_CVISI_baseline = cell(1,3);
    Rate_CVISI_Voc = cell(1,3);
    Rate_CVISI_Trill = cell(1,3);

    if strcmp(SessionName, 'Operant')
        figure(1)
    else
        figure(4)
    end
    clf
    subplot(1,3,1)
    H=histogram2(LogRate_baseline, CV_ISI_baseline, 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5],'NumBins',[24 36]);
    Rate_CVISI_baseline{1} = H.Values;
    Rate_CVISI_baseline{2} = H.XBinEdges;
    Rate_CVISI_baseline{3} = H.YBinEdges;
    set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
    xlabel('Spike Rate (Hz)')
    ylabel('CV of ISI')
    title('Time points before vocalizations (baseline)')

    subplot(1,3,2)
    H=histogram2(LogRate_Voc, CV_ISI_Voc, 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5], 'XBinEdges',H.XBinEdges, 'YBinEdges',H.YBinEdges);
    Rate_CVISI_Voc{1} = H.Values;
    Rate_CVISI_Voc{2} = H.XBinEdges;
    Rate_CVISI_Voc{3} = H.YBinEdges;
    set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
    xlabel('Spike Rate (Hz)')
    ylabel('CV of ISI')
    title('Time points around vocalizations')

    subplot(1,3,3)
    H=histogram2(LogRate_Voc(logical(Trill01)), CV_ISI_Voc(logical(Trill01)), 'FaceColor','flat', 'XBinLimits', [-1 2],'YBinLimits', [0 5], 'XBinEdges',H.XBinEdges, 'YBinEdges',H.YBinEdges);
    Rate_CVISI_Trill{1} = H.Values;
    Rate_CVISI_Trill{2} = H.XBinEdges;
    Rate_CVISI_Trill{3} = H.YBinEdges;
    set(gca,'XTick',log10([0.1 1:1:10 20:10:100]), 'XTickLabel',[0:1:10 20:10:100])
    xlabel('Spike Rate (Hz)')
    ylabel('CV of ISI')
    title('Time points around Trills')
    suplabel(sprintf('%d/%d %s', ss, NCells,SessionName));
    
    
    % if strcmp(SessionName, 'Operant')
    %     figure(2)
    % else
    %     figure(5)
    % end
    % clf
    
    Spectrum_Win_Voc_out = cell(2,1);
    Spectrum_Win_Trill_out = cell(2,1);
    Spectrum_Win_NonTrill_out = cell(2,1);
    Spectrum_Win_baseline_out = cell(2,1);
    if ~isempty(Spectrum_Win_Voc)
        Spectrum_Win_Voc_out{2} = Spectrum_Win_Voc{1,2};
        Spectrum_Win_Voc_out{1} = mean(cell2mat(Spectrum_Win_Voc(:,1)));
        if sum(Trill01)
            Spectrum_Win_Trill_out{2} = Spectrum_Win_Voc{1,2};
            Spectrum_Win_Trill_out{1} = mean(cell2mat(Spectrum_Win_Voc(logical(Trill01),1)));
            Spectrum_Win_NonTrill_out{2} = Spectrum_Win_Voc{1,2};
            Spectrum_Win_NonTrill_out{1} = mean(cell2mat(Spectrum_Win_Voc(logical(~Trill01),1)));
        end
    end
    if ~isempty(Spectrum_Win_Baseline)
        Spectrum_Win_baseline_out{2} = Spectrum_Win_Baseline{1,2};
        Spectrum_Win_baseline_out{1} = mean(cell2mat(Spectrum_Win_Baseline(:,1)));
    end

    % clf
    % if ~isempty(Spectrum_Win_Voc_out)
    %     plot(Spectrum_Win_Voc_out{2}, 10*log10(Spectrum_Win_Voc_out{1}),  '-','Color',[0.8500 0.3250 0.0980], 'LineWidth',2)
    %     hold on
    %     if sum(Trill01)
    %         plot(Spectrum_Win_Trill_out{2}, 10*log10(Spectrum_Win_Trill_out{1}),  '-','Color',[0.9290 0.6940 0.1250], 'LineWidth',2)
    %         hold on
    %         plot(Spectrum_Win_NonTrill_out{2}, 10*log10(Spectrum_Win_NonTrill_out{1}), '-', 'Color', [0.4940 0.1840 0.5560], 'LineWidth',2)
    %         hold on
    %     end
    % end
    % if ~isempty(Spectrum_Win_baseline_out)
    %     plot(Spectrum_Win_baseline_out{2}, 10*log10(Spectrum_Win_baseline_out{1}), '-', 'Color', [0 0.4470 0.7410], 'LineWidth',2)
    % end
    % xlim([0 30]);
    % hold off
    % ylabel('Power of the Spike pattern (DB/Hz)')
    % xlabel('Frequency (Hz)')
    % title(sprintf('%d/%d Spike Power Spectrum %s', ss, NCells, SessionName))
    % if ~isempty(Spectrum_Win_Voc_out) && ~isempty(Spectrum_Win_baseline_out) && any(Trill01)
    %     legend('all Voc','Trills','NonTrills', 'baseline')
    % elseif isempty(Spectrum_Win_Voc_out) && ~isempty(Spectrum_Win_baseline_out)
    %     legend('baseline')
    % elseif ~isempty(Spectrum_Win_Voc_out) && ~isempty(Spectrum_Win_baseline_out) && ~any(Trill01)
    %     legend('all Voc', 'baseline')
    % end

    if strcmp(SessionName, 'Operant')
        figure(3)
    else
        figure(6)
    end
    clf
    histogram([ISI_ms_Behav{~Trill01}], 0:5:max([ISI_ms_Behav{~Trill01}]), 'normalization', 'probability')
    hold on
    histogram([ISI_ms_Behav{logical(Trill01)}], 0:5:max([ISI_ms_Behav{logical(Trill01)}]), 'normalization', 'probability')
    hold on
    histogram([ISI_ms_baseline{logical(Trill01)}], 0:5:max([ISI_ms_baseline{logical(Trill01)}]), 'normalization', 'probability')
    hold on
    histogram([ISI_ms_baseline{~Trill01}], 0:5:max([ISI_ms_baseline{~Trill01}]), 'normalization', 'probability')
    legend('NonTrill','Trill','Baseline Trill', 'Baseline Non Trill')
    ISI_ms.Trill = [ISI_ms_Behav{logical(Trill01)}];
    ISI_ms.NonTrill = [ISI_ms_Behav{~Trill01}];
end


function [Spectrum, Freqs] = spike_pattern_spectrum(Spike_arrival_times_ms, TimeOn, TimeOff)
if nargin<2
    TimeOn = min(Spike_arrival_times_ms);
end
if nargin<3
    TimeOff= max(Spike_arrival_times_ms);
end

% Construct a time varying vector of the neural activity and calculate
    % its spectrum
    % Spike_arrival_times_ms0 = (1 + (Spike_arrival_times - min(Spike_arrival_times)).*10^-3);
    Spike_arrival_times_ms0 = 1 + (Spike_arrival_times_ms - TimeOn);
%     RecordingDuration_s = floor((max(Spike_arrival_times)-min(Spike_arrival_times))/(60*10^6))*60;
%     [KDE,t,Error]=kde_wrapper(Spike_arrival_times_s0,0:0.001:RecordingDuration_s,1000);
    SpikePattern = zeros(round(TimeOff-TimeOn)+1,1);
    for Sp=1:length(Spike_arrival_times_ms0)
        SpikePattern(floor(Spike_arrival_times_ms0(Sp))) = SpikePattern(floor(Spike_arrival_times_ms0(Sp))) +1;
    end
%     SpikePattern = conv(SpikePattern, Expwav,'same');
    if length(SpikePattern)>1024
        [Spectrum, Freqs] = pwelch(SpikePattern - mean(SpikePattern), 1024,512, 1024,1000);
    else
        [Spectrum, Freqs] = pwelch(SpikePattern - mean(SpikePattern), 512,256, 512,1000);
    end
    Spectrum = Spectrum';
    Freqs = Freqs';
end