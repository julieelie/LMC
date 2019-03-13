% Path 2 Data
Path2AllData = '/Volumes/server_home/users/JulieE/LMC_HoHa/logger/';

% Let's find the date of interest in the Google sheet
Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';
[~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:AF200','basic');
Header = RecTableData(1,:);
ND_IDCol = find(contains(Header, 'Neural data voc in operant + PSTH'));
Dates_Ind = find(RecTableData(:,ND_IDCol)==2);
Dates = RecTableData(Dates_Ind,1);

%% Loop through dates and gather data
% For each day, each logger, each unit find the average spike rate
NDays = length(Dates);
HypNSU = 100;
TimeStep = 10*60; % Time resolution in seconds at which the spike rate shoule be calculated
TimePoints = 0:TimeStep:60*60*7;
SpikeRate = nan(HypNSU,length(TimePoints)-1);% will contain spike rate in Hz
Peak2Peak = nan(HypNSU,length(TimePoints)-1);% will contain mean spike amplitude along time
SD_SR = nan(HypNSU,1);
Mean_SR = nan(HypNSU,1);
Bat_ID = nan(HypNSU,1);
Mean_Peak2Peak = nan(HypNSU,1);
SpikeShape = cell(HypNSU,4);
CellCount=0;
for dd=1:NDays
    Date = Dates(dd);
    Path2Data = fullfile(Path2AllData, ['20' Date]);
    % Get the number of loggers
    Logger_dirs = dir(fullfile(Path2Data, '*ogger*'));
    Logger_dirs=Logger_dirs([Logger_dirs.isdir]);
    NLogger = length(Logger_dirs);
    % Identify the type of logger and extract neural data
    LoggerType = cell(NLogger,1);
    for ll=1:NLogger
        LData_folder = fullfile(Logger_dirs(ll).folder, Logger_dirs(ll).name,'extracted_data');
        LDir = dir(fullfile(LData_folder, '*CSC*.mat'));
        LData = load(fullfile(LDir(1).folder, LDir(1).name), 'logger_type', 'logger_serial_number');
        LoggerType{ll}  = LData.logger_type;
        if strcmp(LoggerType{ll}, 'Mous') || strcmp(LoggerType{ll}, 'Rat')
            SU_files = dir(fullfile(LData_folder, '*TT*SS*.mat'));
            NSU = length(SU_files);
            %!!! Get the Bat-ID for that logger!!
            
            % LOAD the data of KDE to do a test for vocalization production
            
            for uu=1:NSU
                CellCount = CellCount+1;
                Cell = load(fullfile(SU_files(uu).folder, SU_files(uu).file)); % This load the Spike_arrival_times in micro seconds and the Spike_snippets in microVolts, ceiled at 500uV for spike sorting projection
                SpikeTimes = (Cell.Spike_arrival_times-Cell.Spike_arrival_times(1))*10^-6; % Spike arrival times centered to the first spike, in seconds
                for tt=1:(length(TimePoints)-1)
                    Spike_local = (SpikeTimes>=TimePoints(tt)) .* (SpikeTimes<=TimePoints(tt+1));
                    SpikeRate(CellCount,tt) = sum(Spike_local)/TimeStep; % this is the local spike rate
                    SpikeInd = find(Spike_local);
                    P2P_local = nan(1,4);
                    for cc=1:4 % calculate for each channel the average peak2peak for the spikes occuring during that time slot
                        P2P_local(cc)=mean(max(Cell.Spike_snippets(SpikeInd,cc,:),[],2) - min(Cell.Spike_snippets(SpikeInd,cc,:),[],2));
                    end
                    Peak2Peak(CellCount,tt) = max(P2P_local);
                end
                % calculate average statistics for that cell
                Mean_Peak2Peak(CellCount) = mean(Peak2Peak(CellCount,:));
                Mean_SR(CellCount) = mean(SpikeRate(CellCount,:));
                SD_SR(CellCount) = std(SpikeRate(CellCount,:));
                Bat_ID(CellCount) = BAT; %% GET BAT ID
                % extract the average spike shape for that cell
                for cc=1:4
                    SpikeShape{CellCount,cc}=nan(2,dim(Cell.Spike_snippets,3));
                    SpikeShape{CellCount,cc}(1,:) = mean(Cell.Spike_snippets(:,cc,:));
                    SpikeShape{CellCount,cc}(2,:) = std(Cell.Spike_snippets(:,cc,:));
                end
                % Test if the cell has a modulated spike rate compared to
                % averaged before vocalizations (one time window with KDE
                % away from the mean by twice SD)
                
            end
            
        end
    end
end

