addpath(genpath('/Users/elie/Documents/CODE/operant_bats'))
addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'))
addpath(genpath('/Users/elie/Documents/CODE/LMC'))
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'))
addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))
Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';

ListOfPaths = {
    %'/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190130/HoHa_190130_1007_VocTrigger_param.txt';
    %'/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190129/HoHa_190129_1023_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190120/HoHa_190120_1208_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190124/HoHa_190124_0957_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190202/HoHa_190202_1046_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190205/HoHa_190205_1140_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190206/HoHa_190206_1024_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190207/HoHa_190207_1136_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190208/HoHa_190208_1018_VocTrigger_param.txt';
    % '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190211/HoHa_190211_1152_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190212/HoHa_190212_1033_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190213/HoHa_190213_1101_VocTrigger_param.txt';
    '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190214/HoHa_190214_1130_VocTrigger_param.txt';
    %'/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190119/HoHa_190119_1158_VocTrigger_param.txt';
    %'/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190118/HoHa_190118_1027_VocTrigger_param.txt';
    %'/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190117/HoHa_190117_1008_VocTrigger_param.txt';
    %'/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190116/HoHa_190116_1126_VocTrigger_param.txt'
    };
for pp=1:length(ListOfPaths)
    Path2ParamFile = ListOfPaths{pp};
    
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190202/HoHa_190202_1046_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190131/HoHa_190131_1108_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190130/HoHa_190130_1007_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190129/HoHa_190129_1023_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190120/HoHa_190120_1208_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190119/HoHa_190119_1158_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190118/HoHa_190118_1027_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190117/HoHa_190117_1008_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190116/HoHa_190116_1126_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190122/HoHa_190122_1027_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190123/HoHa_190123_0943_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190124/HoHa_190124_0957_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190125/HoHa_190125_0925_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190128/HoHa_190128_1022_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190203/HoHa_190203_1259_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190204/HoHa_190204_1051_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190205/HoHa_190205_1140_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190206/HoHa_190206_1024_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190207/HoHa_190207_1136_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190208/HoHa_190208_1018_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190211/HoHa_190211_1152_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190212/HoHa_190212_1033_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190213/HoHa_190213_1101_VocTrigger_param.txt';
%     Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190214/HoHa_190214_1130_VocTrigger_param.txt';
    
    
    
    %% RUN Behavioral and audio data EXtraction
     result_operant_bat(Path2ParamFile)
    
end


%% Loop through neurons to extract spike data
    % Generate the list of paths to gather the data
    BasePath = '/Volumes/server_home/users/JulieE/';
    [ListSSU] = gather_neural_datapath(BasePath);
    
    % Extract the neural data corresponding to the vocalizations identified
    % by voc_localize and voc_localize_operant for each cell
    for ss=1:length(ListSSU)
        cut_neuralData_voc_perfile(InputDataFile, Flags, OutputPath, DenoiseT, Rthreshold)
    end
    
    % Get the path to audio data for operant conditioning experiment
    [AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
    Date = DataFile(6:11);
    fprintf(1,'Working on %s\n', Date)
    ExpStartTime = DataFile(13:16);
    % Set the path to logger data
    Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'logger',['20' Date]);
    
    % Set the time buffer before vocalizations onset
    BufferBeforeOnset = 200; %ms
    fprintf(' EXTRACTING NEURAL DATA CORRESPONDING TO VOCALIZATIONS \n')
    FlagsExtr = [0 0 1 1]; % FlagsExtr(1)= Raw data, FlagsExtr(2) = LFP, FlagsExtr(3) = Tetrodes, FlagsExtr(4) = single units
    DenoiseT=1;
    Rthreshold = [0.92 0.94 0.96 0.98];
    cut_neuralData_voc(Logger_dir,Date, ExpStartTime,FlagsExtr,BufferBeforeOnset,DenoiseT,Rthreshold);
    
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
        Flags=[0 1];% Flags = whether to calculate PSTH of Tetrode (Flags(1)=1) and/or Single units
        % (Flags(2)=1))
        KDE_Cal = 1;
        PLOT=1; % set to 1 to plot the results, 0 to just return data
        fprintf(1,' PSTH of NEURAL DATA CORRESPONDING TO VOCALIZATIONS DURING OPERANT %s\n', NeuroLoggerID)
%         for rr=1:length(Rthreshold)
%             [SpikeTimesVoc.(NeuroLoggerID)]=plot_psth_voc(Logger_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, BufferBeforeOnset, KDE_Cal,PLOT,rr);
%         end
        [SpikeTimesVoc.(NeuroLoggerID)]=plot_psth_voc(Logger_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, BufferBeforeOnset, KDE_Cal,PLOT);
        [SpikeTimesVocOff.(NeuroLoggerID)]=plot_psth_voc_off(Logger_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, BufferBeforeOnset, KDE_Cal,PLOT);
    end
    save(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, BufferBeforeOnset)), 'SpikeTimesVoc','-append')
    save(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, BufferBeforeOnset)), 'SpikeTimesVocOff','-append')
    close all
    %     pause()
    
     %% Extract data of the bats doing other actions during the free behavior session (RecOnly)
    fprintf(' EXTRACTING ONSET/OFFSET TIMES OF OTHER BEHAVIORS DURING FREE SESSION \n')
    RecOnlySession = dir(fullfile(AudioDataPath, '*RecOnly_events.txt'));
    if isempty(RecOnlySession)
        fprintf(1,'No free interaction session on that day!\n')
    else
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
        %
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
        %
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
        MaxDur = 700; % duration by which each long behavioral sequence should be cut often set at 600 but for 190118, 700 is better
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
            PLOT = 0; % set to 1 to plot the results, 0 to just return data
            KDE_Cal = 1;
            fprintf(1,' PSTH of NEURAL DATA CORRESPONDING TO OTHER BEHAVIORS DURING FREE SOCIALIZATION %s \n', NeuroLoggerID)
            [SpikeTimesBehav.(NeuroLoggerID)]= plot_psth_behav(Logger_dir, Date, ExpStartTime, NeuroLoggerID,Bat_ID, Flags, MaxDur, KDE_Cal, PLOT);
        end
        save(fullfile(Logger_dir, sprintf('%s_%s_BehavExtractData_%d.mat', Date, ExpStartTime,BufferBeforeBehavOnset)),'SpikeTimesBehav','-append');
        close all
    end
%     
    %% Plot one PSTH per unit with all actions
    if isempty(RecOnlySession)
            fprintf(1,'No free interaction session on that day!\n')
            for nl=1:length(NL_ID)
                NeuroLoggerID = ['Logger' num2str(NL_ID(nl))];
                SpikeTimesBehav.(NeuroLoggerID) = [];
            end
    else
        load(fullfile(Logger_dir, sprintf('%s_%s_BehavExtractData_%d.mat', Date, ExpStartTime,BufferBeforeBehavOnset)) , 'SpikeTimesBehav');
    end
    for nl=1:length(NL_ID)
        NeuroLoggerID = ['Logger' num2str(NL_ID(nl))];
        Bat_ID = DataInfo{BatIDCol(find(BatIDCol<NLCol(nl),1,'last'))};
        Flags=[1 1];% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
        % (Flags(2)=1))
        KDE_Cal = 1;
        
        fprintf(1,' PSTH of NEURAL DATA CORRESPONDING TO BEHAVIORS DURING FREE SOCIALIZATION AND VOCAL ACTIVITY DURING OPERANT CONDITIONING %s\n', NeuroLoggerID)
        plot_psth_voc_and_behav(SpikeTimesBehav.(NeuroLoggerID),SpikeTimesVoc.(NeuroLoggerID),Logger_dir,Date, NeuroLoggerID,Flags, BufferBeforeOnset,MaxDur, KDE_Cal);
        plot_psth_voc_and_behav_off(SpikeTimesBehav.(NeuroLoggerID),SpikeTimesVocOff.(NeuroLoggerID),Logger_dir,Date, NeuroLoggerID,Flags, BufferBeforeOnset,MaxDur, KDE_Cal);
    end

    
%%   INTERNAL FUNCTIONS 
function [ListSSU] = gather_neural_datapath(BasePath)
fprintf(1,'*** Gathering paths to spike sorted units ***')
ListSSU = cell(10^3,1); % initialize the list to 1000
ExpFolders = dir(fullfile(BasePath,'LMC*'));
NSSU = 0; % counter for single units
for ee=1:length(ExpFolders)
    fprintf(1, '\n  -> Looking into  %s... ', fullfile(ExpFolders(ee).folder,ExpFolders(ee).name))
    DateFolders = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name, 'logger','20*'));
    for dd=1:length(DateFolders)
        fprintf(1, '%s   ', DateFolders(dd).name);
        LoggerFolders = dir(fullfile(DateFolders(dd).folder, DateFolders(dd).name,'Logger*'));
        for ll = 1:length(LoggerFolders)
            SSFiles = dir(LoggerFolders(ll).folder, LoggerFolders(ll).name, 'extracted_data', '*_TT*_SS*.mat');
            if isempty(SSFiles)
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