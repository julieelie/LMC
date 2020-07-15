addpath(genpath('/Users/elie/Documents/CODE/operant_bats'))
addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'))
addpath(genpath('/Users/elie/Documents/CODE/LMC'))
addpath(genpath('/Users/elie/Documents/CODE/Kilosort2'))
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'))
addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))
Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';
BasePath = '/Volumes/Julie4T';

%% Generate the list of paths to gather the raw data
% BasePath = '/Volumes/server_home/users/JulieE/LMC';
OutputPath = '/Users/elie/Documents/LMCResults';
[ListSSU] = gather_neural_datapath(BasePath);
%% Listing datacells that are extracted
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
%% Explore the regimes of each cell
NeuralWin = 1000; %duration of the time window in
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
Buffer = 2000;
NCells = length(CellsPath);
Rate_CVISI_All = cell(NCells,3);
Rate_CVISI_Voc = cell(NCells,3);
Rate_CVISI_NonVoc = cell(NCells,3);
Rate_CVISI_Trill = cell(NCells,3);
Rate_CVISI_Ba = cell(NCells,3);

for ss=1:length(CellsPath)
    fprintf(1,'File %d/%d: %s \n',ss,length(CellsPath),CellsPath(ss).name)
    load(fullfile(CellsPath(ss).folder, CellsPath(ss).name),'FreeBehavSession','OperantSession','PlayBackSession','QualitySSU','What','VocRank')
    Ind_ = strfind(CellsPath(ss).name, '_');
    SUBJ = CellsPath(ss).name(1:(Ind_(1)-1));
    Date = CellsPath(ss).name(Ind_(1)+ (1:8));
    Tetrode = CellsPath(ss).name(Ind_(3)+ 1);
    SSQ = CellsPath(ss).name(Ind_(3) - 1);
    CellID = CellsPath(ss).name((strfind(CellsPath(ss).name, '-') + 1) : (strfind(CellsPath(ss).name, '.') - 1));
    ssu = find(contains(ListSSU, sprintf('%s_%s_TT%s_SS%s_%s',SUBJ, Date, Tetrode, SSQ, CellID)));
    load(ListSSU{ssu}, 'Spike_arrival_times'); % load spike arrival times in useconds
    
    % Find vocalizations onset/offsets
    [Path2Data, DataFile]=fileparts(ListSSU{ssu});
    PathParts = regexp(Path2Data, '/', 'split');
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
    ISI_ms = diff(Spike_arrival_times_ms);
    TimeOn = -2*NeuralWin:Overlap:(max(Spike_arrival_times_ms)-NeuralWin);
    TimeOff = -1*NeuralWin:Overlap:max(Spike_arrival_times_ms);
    Rate = nan(length(TimeOn),1);
    STD_ISI = nan(length(TimeOn),1);
    Mean_ISI = nan(length(TimeOn),1);
    Voc01 = nan(length(TimeOn),1);
    if sum(contains(VocRank, 'first')) == size(VocDataTime,1) % we can retrieve the type of call!
        Trill01 = zeros(length(TimeOn),1); % Trill =1; Non-Trill or non Voc = 0 -> Non-Trill = Voc01 .* ~Trill01
    end
    for bb=1:length(TimeOn)
        Spikes01 = (Spike_arrival_times_ms>=TimeOn(bb)) .* (Spike_arrival_times_ms<TimeOff(bb));
        Rate(bb) = sum(Spikes01)/(NeuralWin*10^-3);
        STD_ISI(bb) = std(ISI_ms(logical(Spikes01(1:(end-1)))));
        Mean_ISI(bb) = mean(ISI_ms(logical(Spikes01(1:(end-1)))));
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
    CV_ISI = STD_ISI ./ Mean_ISI;
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
    Voc01 = Voc01(1:(FirstDead-2));
    if sum(contains(VocRank, 'first')) == size(VocDataTime,1)
        Trill01 = Trill01(1:(FirstDead-2));
    end
    
    
    % Construct a time varying vector of the neural activity and calculate
    % its spectrum
    Spike_arrival_times_ms0 = (1 + (Spike_arrival_times - min(Spike_arrival_times)).*10^-3);
%     RecordingDuration_s = floor((max(Spike_arrival_times)-min(Spike_arrival_times))/(60*10^6))*60;
%     [KDE,t,Error]=kde_wrapper(Spike_arrival_times_s0,0:0.001:RecordingDuration_s,1000);
    
    SpikePattern = zeros(round(max(Spike_arrival_times_ms0)),1);
    for Sp=1:length(Spike_arrival_times_ms0)
        SpikePattern(floor(Spike_arrival_times_ms0(Sp))) = SpikePattern(floor(Spike_arrival_times_ms0(Sp))) +1;
    end
%     SpikePattern = conv(SpikePattern, Expwav,'same');
    [Spectrum, Freqs] = pwelch(SpikePattern - mean(SpikePattern), 4096,2048, 4096,1000);
    
    
    if sum(contains(VocRank, 'first')) == size(VocDataTime,1)
        figure(1)
        clf
        H=histogram2(log10(Rate), CV_ISI, 'FaceColor','flat', 'XBinLimits', [0 2],'YBinLimits', [0 5], 'NumBins',[24 36]);
        XBinEdges = H.XBinEdges;
        YBinEdges = H.YBinEdges;
        Rate_CVISI_All{ss,1} = H.Values;
        Rate_CVISI_All{ss,2} = H.XBinEdges;
        Rate_CVISI_All{ss,3} = H.YBinEdges;
%         set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
%         xlabel('Spike Rate (Hz)')
%         ylabel('CV of ISI')
%         title('all time points')

        subplot(2,2,1)
        H=histogram2(log10(Rate(logical(Voc01))), CV_ISI(logical(Voc01)), 'FaceColor','flat', 'XBinLimits', [0 2],'YBinLimits', [0 5], 'XBinEdges',XBinEdges, 'YBinEdges',YBinEdges);
        Rate_CVISI_Voc{ss,1} = H.Values;
        Rate_CVISI_Voc{ss,2} = H.XBinEdges;
        Rate_CVISI_Voc{ss,3} = H.YBinEdges;
        set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title(sprintf('Time points around vocalizations (within %d ms)', Buffer+Overlap))

        subplot(2,2,2)
        H=histogram2(log10(Rate(~Voc01)), CV_ISI(~Voc01), 'FaceColor','flat', 'XBinLimits', [0 2],'YBinLimits', [0 5], 'XBinEdges',XBinEdges, 'YBinEdges',YBinEdges);
        Rate_CVISI_NonVoc{ss,1} = H.Values;
        Rate_CVISI_NonVoc{ss,2} = H.XBinEdges;
        Rate_CVISI_NonVoc{ss,3} = H.YBinEdges;
        set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title('Time points outside of vocalizations')
        
        subplot(2,2,3)
        H=histogram2(log10(Rate(logical(Voc01.*Trill01))), CV_ISI(logical(Voc01.*Trill01)), 'FaceColor','flat', 'XBinLimits', [0 2],'YBinLimits', [0 5], 'XBinEdges',XBinEdges, 'YBinEdges',YBinEdges);
        Rate_CVISI_Trill{ss,1} = H.Values;
        Rate_CVISI_Trill{ss,2} = H.XBinEdges;
        Rate_CVISI_Trill{ss,3} = H.YBinEdges;
        set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title(sprintf('Time points around Trills (within %d ms)', Buffer+Overlap))
        
        subplot(2,2,4)
        H=histogram2(log10(Rate(logical(Voc01.*~Trill01))), CV_ISI(logical(Voc01.*~Trill01)), 'FaceColor','flat', 'XBinLimits', [0 2],'YBinLimits', [0 5], 'XBinEdges',XBinEdges, 'YBinEdges',YBinEdges);
        Rate_CVISI_Ba{ss,1} = H.Values;
        Rate_CVISI_Ba{ss,2} = H.XBinEdges;
        Rate_CVISI_Ba{ss,3} = H.YBinEdges;
        set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title(sprintf('Time points around Barks (within %d ms)', Buffer+Overlap))
    else
        figure(1)
        clf
        subplot(1,3,1)
        H=histogram2(log10(Rate), CV_ISI, 'FaceColor','flat', 'XBinLimits', [0 2],'YBinLimits', [0 5],'NumBins',[24 36]);
        Rate_CVISI_All{ss,1} = H.Values;
        Rate_CVISI_All{ss,2} = H.XBinEdges;
        Rate_CVISI_All{ss,3} = H.YBinEdges;
        set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title('all time points')

        subplot(1,3,2)
        H=histogram2(log10(Rate(logical(Voc01))), CV_ISI(logical(Voc01)), 'FaceColor','flat', 'XBinLimits', [0 2],'YBinLimits', [0 5], 'XBinEdges',H.XBinEdges, 'YBinEdges',H.YBinEdges);
        Rate_CVISI_Voc{ss,1} = H.Values;
        Rate_CVISI_Voc{ss,2} = H.XBinEdges;
        Rate_CVISI_Voc{ss,3} = H.YBinEdges;
        set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title(sprintf('Time points around vocalizations (within %d ms)', Buffer+Overlap))

        subplot(1,3,3)
        H=histogram2(log10(Rate(~Voc01)), CV_ISI(~Voc01), 'FaceColor','flat', 'XBinLimits', [0 2],'YBinLimits', [0 5], 'XBinEdges',H.XBinEdges, 'YBinEdges',H.YBinEdges);
        Rate_CVISI_NonVoc{ss,1} = H.Values;
        Rate_CVISI_NonVoc{ss,2} = H.XBinEdges;
        Rate_CVISI_NonVoc{ss,3} = H.YBinEdges;
        set(gca,'XTick',log10([1:1:10 20:10:100]), 'XTickLabel',[1:1:10 20:10:100])
        xlabel('Spike Rate (Hz)')
        ylabel('CV of ISI')
        title('Time points outside of vocalizations')
    end
    
    
    
    figure(2)
    clf
    subplot(2,1,1)
    imagesc((TimeOn(1:FirstDead-1)+NeuralWin/2).*10^-3,0,Rate')
    xlabel('Time (s)')
    zlabel('Hz')
    title('Spike Rate (Hz)')
    colorbar()
    subplot(2,1,2)
    imagesc((TimeOn(1:FirstDead-1)+NeuralWin/2).*10^-3,0,CV_ISI')
    xlabel('Time (s)')
    title('CV ISI')
    colorbar()
    
    figure(3)
    clf
    plot(Freqs,log(Spectrum), '-k', 'LineWidth',2)
    xlim([0 100]);
    ylabel('Log Power of the Spike Rate')
    xlabel('Frequency (Hz)')
    pause(1)
end
fprintf(' DONE \n')
save(fullfile(OutputPath, 'CellSpikingRegime.mat'),'ListSSU', 'CellsPath', 'NeuralWin','Overlap', 'Rate_CVISI_All', 'Rate_CVISI_Voc', 'Rate_CVISI_NonVoc', 'Rate_CVISI_Trill', 'Rate_CVISI_Ba');


%% Run a umap on the CV_Rate plots
RateCVISI_mat = nan(size(Rate_CVISI_All,1), numel(Rate_CVISI_All{1,1}));
for cc=1:NCells
    RateCVISI_mat(cc,:) = reshape(Rate_CVISI_All{cc,1}, 1, numel(Rate_CVISI_All{cc,1}));
end
figure()
imagesc(RateCVISI_mat)
xlabel('CV on ISI and Rate profile')
ylabel('Cell#')
[Reduction,UMAP,ClustID]= run_umap(RateCVISI_mat);


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