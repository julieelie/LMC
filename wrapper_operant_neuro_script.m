%Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190202/HoHa_190202_1046_VocTrigger_param.txt';
% Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190131/HoHa_190131_1108_VocTrigger_param.txt';
% Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190130/HoHa_190130_1007_VocTrigger_param.txt';
Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190129/HoHa_190129_1023_VocTrigger_param.txt';
% Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190201/HoHa_190201_1023_VocTrigger_param.txt';

addpath(genpath('/Users/elie/Documents/CODE/operant_bats'))
addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'))
addpath(genpath('/Users/elie/Documents/CODE/LMC'))
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'))
addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))

%% Some important parameters, inputs
% Get the path to audio data for operant conditioning experiment
[AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
Date = DataFile(6:11);
ExpStartTime = DataFile(13:16);

Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';

% Set the path to logger data
Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'logger',['20' Date]);

% Set the time buffer before vocalizations onset
BufferBeforeOnset = 200; %ms

%% RUN Behavioral and audio data EXtraction
result_operant_bat(Path2ParamFile)

%% Extract the neural data corresponding to the vocalizations
fprintf(' EXTRACTING NEURAL DATA CORRESPONDING TO VOCALIZATIONS \n')
FlagsExtr = [0 0 1 1]; % FlagsExtr(1)= Raw data, FlagsExtr(2) = LFP, FlagsExtr(3) = Tetrodes, FlagsExtr(4) = single units
cut_neuralData_voc(Logger_dir,Date, ExpStartTime,FlagsExtr,BufferBeforeOnset);

%% Plot PSTH of the bats hearing or producing a vocalization during the operant conditioning
% Find the ID of the Neural loggers and corresponding audiologger for each implanted bat
[~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:P200','basic');
RowData = find((cell2mat(RecTableData(2:end,1))== str2double(Date))) +1;
DataInfo = RecTableData(RowData,:);
Header = RecTableData(1,:);
BatIDCol = find(contains(Header, 'Bat'));
NLCol = find(contains(Header, 'NL'));
ALCol = find(contains(Header, 'AL-throat'));
NL_ID = cell2mat(DataInfo(NLCol));
for nl=1:length(NL_ID)
    NeuroLoggerID = ['Logger' num2str(NL_ID(nl))];
    AL_ID = DataInfo{ALCol(find(ALCol<NLCol(nl),1,'last'))};
    AudioLoggerID = ['Logger' num2str(AL_ID)];
    Flags=[1 1];% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
% (Flags(2)=1))
    KDE_Cal = 1;
    fprintf(' PSTH of NEURAL DATA CORRESPONDING TO VOCALIZATIONS DURING OPERANT \n')
    [SpikeTimesVoc.(NeuroLoggerID)]=plot_psth_voc(Logger_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, BufferBeforeOnset, KDE_Cal);
end
close all
%     pause()

%% Extract data of the bats doing other actions during the free behavior session (RecOnly)
fprintf(' EXTRACTING ONSET/OFFSET TIMES OF OTHER BEHAVIORS DURING FREE SESSION \n')
RecOnlySession = dir(fullfile(AudioDataPath, '*RecOnly_events.txt'));
if length(RecOnlySession)>1
    fprintf(1, 'Several RecOnly session were done on that day:\n')
    for ss=1:length(RecOnlySession)
        fprintf(1, '%d. %s\n', ss, RecOnlySession(ss).name);
    end
    Inputss = input('Your choice:\n');
    RecOnlySession = RecOnlySession(Inputss);
end
Date = RecOnlySession.name(6:11);
ExpStartTime = RecOnlySession.name(13:16);
% extract the time onset/offset of behaviors
get_logger_data_behav(AudioDataPath, Logger_dir, Date, ExpStartTime) 

%% Extract the neural data corresponding to other actions during the free behavior session (RecOnly)
fprintf(' EXTRACTING NEURAL DATA CORRESPONDING TO OTHER BEHAVIORS DURING FREE SESSION \n')
FlagsExtr = [0 0 1 1]; % FlagsExtr(1)= Raw data, FlagsExtr(2) = LFP, FlagsExtr(3) = Tetrodes, FlagsExtr(4) = single units
BufferBeforeBehavOnset = 0;
RecOnlySession = dir(fullfile(AudioDataPath, '*RecOnly_events.txt'));
if length(RecOnlySession)>1
    fprintf(1, 'Several RecOnly session were done on that day:\n')
    for ss=1:length(RecOnlySession)
        fprintf(1, '%d. %s\n', ss, RecOnlySession(ss).name);
    end
    Inputss = input('Your choice:\n');
    RecOnlySession = RecOnlySession(Inputss);
end
Date = RecOnlySession.name(6:11);
ExpStartTime = RecOnlySession.name(13:16);
cut_neuralData_behav(Logger_dir,Date, ExpStartTime,FlagsExtr,BufferBeforeBehavOnset);

 %% Plot PSTH of the bats doing other actions during free socialization!
 RecOnlySession = dir(fullfile(AudioDataPath, '*RecOnly_events.txt'));
if length(RecOnlySession)>1
    fprintf(1, 'Several RecOnly session were done on that day:\n')
    for ss=1:length(RecOnlySession)
        fprintf(1, '%d. %s\n', ss, RecOnlySession(ss).name);
    end
    Inputss = input('Your choice:\n');
    RecOnlySession = RecOnlySession(Inputss);
end
Date = RecOnlySession.name(6:11);
ExpStartTime = RecOnlySession.name(13:16);
 % Find the ID of the Neural loggers and corresponding ID tag of each implanted bat
MaxDur = 600; % duration by which each long behavioral sequence should be cut
[~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:P200','basic');
RowData = find((cell2mat(RecTableData(2:end,1))== str2double(Date))) +1;
DataInfo = RecTableData(RowData,:);
Header = RecTableData(1,:);
BatIDCol = find(contains(Header, 'Bat'));
NLCol = find(contains(Header, 'NL'));
NL_ID = cell2mat(DataInfo(NLCol));
for nl=1:length(NL_ID)
    NeuroLoggerID = ['Logger' num2str(NL_ID(nl))];
    Bat_ID = DataInfo{BatIDCol(find(BatIDCol<NLCol(nl),1,'last'))};
    Flags=[1 1];% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
% (Flags(2)=1))
    KDE_Cal = 1;
    fprintf(' PSTH of NEURAL DATA CORRESPONDING TO OTHER BEHAVIORS DURING FREE SOCIALIZATION \n')
    [SpikeTimesBehav.(NeuroLoggerID)]= plot_psth_behav(Logger_dir, Date, ExpStartTime, NeuroLoggerID,Bat_ID, Flags, MaxDur, KDE_Cal);
end
close all

%% Plot one PSTH per unit with all actions
for nl=1:length(NL_ID)
    NeuroLoggerID = ['Logger' num2str(NL_ID(nl))];
    Bat_ID = DataInfo{BatIDCol(find(BatIDCol<NLCol(nl),1,'last'))};
    Flags=[1 1];% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
% (Flags(2)=1))
    KDE_Cal = 1;
    fprintf(' PSTH of NEURAL DATA CORRESPONDING TO BEHAVIORS DURING FREE SOCIALIZATION AND VOCAL ACTIVITY DURING OPERANT CONDITIONING \n')
    plot_psth_voc_and_behav(SpikeTimesBehav.(NeuroLoggerID),SpikeTimesVoc.(NeuroLoggerID),Logger_dir,Date, NeuroLoggerID,Flags, BufferBeforeOnset, KDE_Cal);
end