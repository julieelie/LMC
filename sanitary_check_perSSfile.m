function  sanitary_check_perSSfile(InputDataFile, OutputPath)

% INPUT:
%       InputDataFile: path to a single unit mat file.

%       OutputPath: string that indicate where the data should be saved. If
%       not specified or left empty, the data will be saved in the folder
%       of the inputfile

%% Hard coded calculation features
TimeStep = 1*60; % Time resolution in seconds at which the spike rate shoule be calculated
RateThreshold = 1/(2*TimeStep); % Frequency in Hertz below which the cell is considered as dead if happening over 5 steps (5min 0 spikes), it conrresponds to default value for the rate in kde_wrapper when there is no spike

%% Massage the input and get info about the SS data
% Decide about the output directory
[Path2Data, DataFile]=fileparts(InputDataFile);
PathParts = regexp(Path2Data, '/', 'split');
Loggers_dir = ['/' fullfile(PathParts{1:end-2})];
if nargin<2 || isempty(OutputPath)
    OutputPath = Loggers_dir;
end

% Get the date of the recording
Idx_ = strfind(DataFile, '_');
Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));

% Get the tetrode ID
NeuralInputID{1} = DataFile(strfind(DataFile, 'TT')+2);
% Get the SS ID
NeuralInputID{2} = DataFile((Idx_(end)+1):end);
% Get the SS quality
NeuralInputID{3} = DataFile(strfind(DataFile, 'SS')+2);

% Get the subject ID
SubjectID = DataFile(1:5);

% Get the corresponding tetrode file
TetrodeFile = dir(fullfile(Path2Data, sprintf('%s_%s_TT%s*Sorted*.mat',SubjectID,Date,NeuralInputID{1})));
if isempty(TetrodeFile)
    TetrodeFile = dir(fullfile(Path2Data, sprintf('%s_%s_TTAll_SSorted.mat',SubjectID,Date)));
end

%% Get the time in transceiver time of each session in s
All_loggers = dir(fullfile(Loggers_dir, '*ogger*'));
DirFlags = [All_loggers.isdir];
% Extract only those that are directories.
All_loggers = All_loggers(DirFlags);
NLog = length(All_loggers);
FreeBehavSession = nan(1,2); % Onset/offset in transceiver time in s of the experiment Rec Only
OperantSession = nan(1,2);% Onset/offset in transceiver time in s of the experiment All Voc Reward (operant conditioning)
PlayBackSession = nan(1,2);% Onset/offset in transceiver time in s of the experiment Playback
for ll=1:NLog
%     fprintf('%s ', All_loggers(ll).name);
    Eventfile = dir(fullfile(All_loggers(ll).folder,All_loggers(ll).name, 'extracted_data', '*_EVENTS.mat')); % load file with TTL status info
    load(fullfile(Eventfile.folder, Eventfile.name), 'event_types_and_details', 'event_timestamps_usec');
    B=find(cellfun(@(x) contains(x,'all voc reward start'),event_types_and_details));
    C=find(cellfun(@(x) contains(x,'all voc reward stop'),event_types_and_details));
    D = find(cellfun(@(x) contains(x,'rec only start'),event_types_and_details));
    E = find(cellfun(@(x) contains(x,'rec only stop'),event_types_and_details));
    F = find(cellfun(@(x) contains(x,'playback start'),event_types_and_details));
    G = find(cellfun(@(x) contains(x,'playback stop'),event_types_and_details));
    if ~isempty(B)
        if length(B)>1
            fprintf(1,'There are %d all-voc-reward-start\n', length(B))
            for ii=1:length(B)
                1e-6*event_timestamps_usec(B(ii))
            end
            Ii = input('Which one do you want to choose?\n');
            B = B(Ii);
        end
        OperantSession(1) = 1e-6*event_timestamps_usec(B);
    end
    if ~isempty(C)
        if length(C)>1
            fprintf(1,'There are %d all-voc-reward-stop\n', length(C))
            for ii=1:length(C)
                1e-6*event_timestamps_usec(C(ii))
            end
            Ii = input('Which one do you want to choose?\n');
            C = C(Ii);
        end
        OperantSession(2) = 1e-6*event_timestamps_usec(C);
    end
    if ~isempty(D)
        if length(D)>1
            fprintf(1,'There are %d rec-only-start\n', length(D))
            for ii=1:length(D)
                1e-6*event_timestamps_usec(D(ii))
            end
            Ii = input('Which one do you want to choose?\n');
            D = D(Ii);
        end
        FreeBehavSession(1) = 1e-6*event_timestamps_usec(D);
    end
    if ~isempty(E)
        if length(E)>1
            fprintf(1,'There are %d rec-only-stop\n', length(E))
            for ii=1:length(E)
                1e-6*event_timestamps_usec(E(ii))
            end
            Ii = input('Which one do you want to choose?\n');
            E = E(Ii);
        end
        FreeBehavSession(2) = 1e-6*event_timestamps_usec(E);
    end
    if ~isempty(F)
        if length(F)>1
            fprintf(1,'There are %d playback-start\n', length(F))
            for ii=1:length(F)
                1e-6*event_timestamps_usec(F(ii))
            end
            Ii = input('Which one do you want to choose?\n');
            F = F(Ii);
        end
        PlayBackSession(1) = 1e-6*event_timestamps_usec(F);
    end
    if ~isempty(G)
        if length(G)>1
            fprintf(1,'There are %d playback-stop\n', length(G))
            for ii=1:length(G)
                1e-6*event_timestamps_usec(G(ii))
            end
            Ii = input('Which one do you want to choose?\n');
            G = G(Ii);
        end
        PlayBackSession(2) = 1e-6*event_timestamps_usec(G);
    end
end

if isnan(OperantSession(2)) && ~isnan(FreeBehavSession(1))
    OperantSession(2) = FreeBehavSession(1);
elseif ~isnan(OperantSession(2)) && isnan(FreeBehavSession(1))
    FreeBehavSession(1) = OperantSession(2);
end
if isnan(FreeBehavSession(2)) && ~isnan(PlayBackSession(1))
    FreeBehavSession(2) = PlayBackSession(1);
elseif ~isnan(FreeBehavSession(2)) && isnan(PlayBackSession(1))
    PlayBackSession(1) = FreeBehavSession(2);
end


%% Load the data and calculate KDE of the rate
Cell = load(fullfile(InputDataFile)); % This load the Spike_arrival_times in micro seconds and the Spike_snippets in microVolts, ceiled at 500uV for spike sorting projection
if isnan(OperantSession(1))
    OperantSession(1) = Cell.Spike_arrival_times(1)*(10^-6);
end
% Define the time points at which KDE should be calculated
MaxRecTime = max([FreeBehavSession OperantSession PlayBackSession])- min([FreeBehavSession OperantSession PlayBackSession]);
MinRecTime = -60*5;
TimePoints = MinRecTime:TimeStep:MaxRecTime;
% Center spikes and calculate KDE
SpikeTimes = (Cell.Spike_arrival_times*(10^-6)-OperantSession(1)); % Spike arrival times (s) centered to the onset of the first experiment or the first spike if the onset of the experiment is unavailable
[KDE,~,KDE_error] = kde_wrapper(SpikeTimes,TimePoints(2:end)-TimeStep/2,1/TimeStep);
% Replace values of spike rate by NaNs for periods of time where the neuron
% was not recorded
% Between Operant and FreeBehav
DeadTime = ((TimePoints(2:end)-TimeStep/2)>diff(OperantSession)) .* ((TimePoints(2:end)-TimeStep/2)<(FreeBehavSession(1)-OperantSession(1)));
KDE(find(DeadTime)) = NaN;
KDE_error(:,find(DeadTime))=NaN;
% Between Free Behav and Play back
if ~isnan(PlayBackSession(2))
    DeadTime = ((TimePoints(2:end)-TimeStep/2)>(FreeBehavSession(2)-OperantSession(1))) .* ((TimePoints(2:end)-TimeStep/2)<(PlayBackSession(1)-OperantSession(1)));
    KDE(find(DeadTime)) = NaN;
    KDE_error(:,find(DeadTime))=NaN;
else
    DeadTime = ((TimePoints(2:end)-TimeStep/2)>(FreeBehavSession(2)-OperantSession(1)));
    KDE(find(DeadTime)) = NaN;
    KDE_error(:,find(DeadTime))=NaN;
end


%% Calculate the peak to peak value for each spike in the best channel or the max
% Find the best channel to calculate the peak2peak value
% (channel with largest spike trace)
P2P_all = nan(size(Cell.Spike_snippets,3),4);
for cc=1:4 % calculate for each channel the peak2peak for all spikes
    P2P_all(:,cc)=reshape(max(Cell.Spike_snippets(:,cc,:),[],1) - min(Cell.Spike_snippets(:,cc,:),[],1),size(Cell.Spike_snippets,3),1);
end
[~,I]=max(P2P_all,[],2);
Best_c = nan(1,4);
for cc=1:4
    Best_c(cc) = sum(I==cc);
end
[~,Best_c] = max(Best_c); % This is the channel with the largest spike
Peak2Peak = squeeze(max(Cell.Spike_snippets(:,Best_c,:),[],1) - min(Cell.Spike_snippets(:,Best_c,:),[],1));
Peak2Peak = Peak2Peak/max(Peak2Peak); % Get a value betwen 0 and 1 for the size of the spike.

%% Plot the KDE along time with the zones for each session
Fig = figure();
if ~isempty(TetrodeFile)
    SS1=subplot(2,1,1);
else
    SS1 = subplot(1,1,1);
end
        
shadedErrorBar((TimePoints(2:end)-TimeStep/2)/60,KDE,KDE_error,{'k-', 'LineWidth',2})
ylabel('Spike rate (Hz)')
xlabel('Time (min)')
MaxYLim = max(SS1.YLim);
hold on
patch(([OperantSession flip(OperantSession)]-OperantSession(1))/60, [0 0 MaxYLim*ones(1,2)],ones(1,4), 'FaceColor', [1 0.8 0], 'FaceAlpha', 0.3);
hold on
fill(([FreeBehavSession flip(FreeBehavSession)]-OperantSession(1))/60, [0 0 MaxYLim*ones(1,2)],ones(1,4), 'FaceColor', [0.2 0.7 1], 'FaceAlpha', 0.3);
hold on
if ~isnan(PlayBackSession(2))
    fill(([PlayBackSession flip(PlayBackSession)]-OperantSession(1))/60, [0 0 MaxYLim*ones(1,2)],ones(1,4), 'FaceColor', [0.3 0.3 0.3], 'FaceAlpha', 0.3);
end
hold on
shadedErrorBar((TimePoints(2:end)-TimeStep/2)/60,KDE,KDE_error,{'k-', 'LineWidth',2})
% Plot the spike arrival times for that cell on top of the previous plot
hold on
% scatter(SpikeTimes/60, MaxYLim/5.*rand(size(SpikeTimes)) + MaxYLim,1,'k')
scatter(SpikeTimes/60, MaxYLim/3.*Peak2Peak + MaxYLim,1,'k')
hold on
title(sprintf('%s on %s TT%s%s SS%s', SubjectID, Date, NeuralInputID{3},NeuralInputID{1},NeuralInputID{2}))

%% get the spike sorting quality for that cell along time

if ~isempty(TetrodeFile)
    TData = load(fullfile(TetrodeFile.folder, TetrodeFile.name));
    SS_i = find(TData.UnitClusters == str2double(NeuralInputID{2}));
    TimeStepQ = diff(TData.TimePoints(1:2));
    XLimSS1 = SS1.XLim;
    SS2=subplot(2,1,2);
    Xtime = ((TData.TimePoints(2:end)-TimeStepQ/2).*10^-6 - OperantSession(1))/60;
    yyaxis left
    plot(Xtime, TData.TimeLRatio(:,SS_i), 'b-', 'LineWidth',2);
    ylabel('LRatio')
    hold on
    yyaxis right
    plot(Xtime, log10(TData.TimeIsolationDistance(:,SS_i)), 'r-', 'LineWidth',2);
    ylabel('IsolationDistance (log10 scale)')
    xlabel('Time (min)')
    SS2.XLim = XLimSS1;
    title(sprintf('cluster %d, LRatio = %.1f, IDist = %.1f', TData.UnitClusters(SS_i), TData.LRatio(SS_i),TData.IsolationDistance(SS_i)));
    hold off
    
end

%% decide of the stability of the cell
DurDead = 10; % dead period are consecutive DurDead minutes below the RateThreshold
% find consecutive periods of 10 min very low rate
DeadCell_Idx = strfind(KDE<=RateThreshold, ones(1,DurDead)); 
if ~isempty(DeadCell_Idx)
    DeadCell_Idx = unique(reshape(repmat(DeadCell_Idx',1,DurDead) + repmat(0:(DurDead-1),length(DeadCell_Idx),1), length(DeadCell_Idx)*DurDead,1));
    % plot these unstable data points
    if ~isempty(TetrodeFile)
        SS1=subplot(2,1,1);
    else
        SS1 = subplot(1,1,1);
    end
    hold on
    plot((TimePoints(DeadCell_Idx+1)-TimeStep/2)/60, KDE(DeadCell_Idx), 'r+')


    DD = find(diff(DeadCell_Idx)>1); % These are the indices that correspond to offset/onset of dead periods
    QualitySSU.DeadTime_usec(1,1) = TimePoints(min(DeadCell_Idx))*60*10^6+ OperantSession(1);
    if ~isempty(DD)
        for dd=1:length(DD)
            QualitySSU.DeadTime_usec(dd,2) = TimePoints(DeadCell_Idx(DD(dd)))*60*10^6+ OperantSession(1);
            QualitySSU.DeadTime_usec(dd,1) = TimePoints(DeadCell_Idx(DD(dd)+1))*60*10^6+ OperantSession(1);
        end
        QualitySSU.DeadTime_usec(dd,2) = TimePoints(min(length(TimePoints), DeadCell_Idx(end)))*60*10^6+ OperantSession(1);
    else
        QualitySSU.DeadTime_usec(1,2) = TimePoints(min(length(TimePoints), DeadCell_Idx(end)))*60*10^6+ OperantSession(1);
    end
else
    QualitySSU.DeadTime_usec = [];
end
%% Stuff in output and save figures and output
QualitySSU.KDE = KDE;
QualitySSU.KDE_error = KDE_error;
QualitySSU.TimePoints = TimePoints;
if ~isempty(TetrodeFile)
    QualitySSU.LRatioIsolationDistTimePoints =TData.TimePoints;
    QualitySSU.TimeLRatio = TData.TimeLRatio(:,SS_i);
    QualitySSU.TimeIsolationDistance = TData.TimeIsolationDistance(:,SS_i);
    QualitySSU.LRatio = TData.LRatio(SS_i);
    QualitySSU.IsolationDistance = TData.TimeIsolationDistance(SS_i);
end

orient(Fig,'landscape')
Fig.PaperPositionMode = 'auto';
set(Fig,'PaperOrientation','landscape');
print(Fig,fullfile(OutputPath,sprintf('%s_%s_SS%s%s-%s_HealthCheck.pdf', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2})),'-dpdf','-fillpage')


OutputFile = fullfile(OutputPath, sprintf('%s_%s_SS%s%s-%s.mat', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2}));
if exist(OutputFile, 'file')
    save(OutputFile, 'QualitySSU','OperantSession','FreeBehavSession','PlayBackSession','-append');
else
    save(OutputFile, 'QualitySSU','OperantSession','FreeBehavSession','PlayBackSession');
end
close all
end