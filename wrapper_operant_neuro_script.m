%Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190202/HoHa_190202_1046_VocTrigger_param.txt';
% Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190131/HoHa_190131_1108_VocTrigger_param.txt';
Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190130/HoHa_190130_1007_VocTrigger_param.txt';
% Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190129/HoHa_190129_1023_VocTrigger_param.txt';

addpath(genpath('/Users/elie/Documents/CODE/operant_bats'))
addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'))
addpath(genpath('/Users/elie/Documents/CODE/LMC'))

%% Some important parameters, inputs
% Get the path to audio data
[AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
Date = DataFile(6:11);
ExpStartTime = DataFile(13:16);

Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';

% Set the path to logger data
Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'logger',['20' Date]);

% Set the time buffer before behavior onset
BufferBeforeOnset = 200; %ms

%% RUN Behavioral and audio data EXtraction
result_operant_bat(Path2ParamFile)

%% Extract the neural data corresponding to the vocalizations
fprintf(' EXTRACTING NEURAL DATA CORRESPONDING TO VOCALIZATIONS \n')
FlagsExtr = [0 1 1 1]; % FlagsExtr(1)= Raw data, FlagsExtr(2) = LFP, FlagsExtr(3) = Tetrodes, FlagsExtr(4) = single units
cut_neuralData_voc(Logger_dir,Date, ExpStartTime,FlagsExtr,BufferBeforeOnset);

%% Plot PSTH of the bats hearing or producing a vocalization
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
    Flags=[1 1];
    KDE_Cal = 0;
    fprintf(' PSTH of NEURAL DATA CORRESPONDING TO VOCALIZATIONS \n')
    plot_psth_voc(Logger_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, BufferBeforeOnset, KDE_Cal)
end
close all
%     pause()

%% Plot PSTH of the bats doing other actions!
addpath(genpath('C:\Users\Julie\Documents\GitHub\GeneralCode\'))
AudioLoggerID = {'Logger5';'Logger5' ; 'Logger7' ;'Logger7';'Logger5';'Logger5';'Logger5';'Logger5';'Logger5';'Logger5';'Logger5'};
NeuroLoggerID = 'Logger16';
Flags=[1 1];
for dd=1:length(Days)
    fprintf(' PSTH of NEURAL DATA CORRESPONDING TO ALL BEHAVIORS \n')
    fprintf('*********** %s *************\n', Days{dd})
    Audio_dir = fullfile(Server_audio_path, Days{dd});
    Loggers_dir = fullfile(Server_logger_path, Days{dd});
    % find the time stamp of each experiment
    StampFiles = dir(fullfile(Audio_dir, '*RecOnly_param.txt'));
    if length(StampFiles)>1
        fprintf('Several Recording Only tests were done on that day, please choose the one you want to look at:\n')
        for ff=1:length(StampFiles)
            fprintf('%d: %s\n',ff, StampFiles(ff).name);
        end
        Indff = input('Your choice:\n');
    else
        Indff = 1;
    end
    ExpStartTime = StampFiles(Indff).name(13:16);
    ExtData_dir = fullfile(Loggers_dir, NeuroLoggerID, 'extracted_data');
    get_logger_data_event(ExtData_dir, Days{dd}, 200, NeuroLoggerID)
    close all
%     pause()
end
