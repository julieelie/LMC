function result_reconly_bat(Path2ParamFile, Path2RecordingTable, Logger_dir)
addpath(genpath('/Users/elie/Documents/CODE/LMC'))
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'))
addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))
TranscExtract = 1; % set to 1 to extract logger data and transceiver time
ForceExtract = 0; % set to 1 to redo the extraction of loggers otherwise the calculations will use the previous extraction data
ForceAllign = 0; % In case the TTL pulses allignment was already done but you want to do it again, set to 1
ForceVocExt1 = 1; % In case the localization on raw files of vocalizations that were manually extracted was already done but you want to do it again set to 1
ForceVocExt2 = 1; % In case the localization on Loggers of vocalizations that were manually extracted was already done but you want to do it again set to 1
ForceWhoID = 1; % In case the identification of bats was already done but you want to re-do it again
close all

% Get the recording date
[AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
Date = DataFile(6:11);

% Set the path to the manual extracts
VocManExtDir = fullfile( '/Users/elie/Google Drive/BatmanData/BatmanCuts', ['20' Date]);

if TranscExtract && nargin<2
    % Set the path to the recording log
    Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';
end
if TranscExtract && nargin<3
    % Set the path to logger data
    Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'logger',['20' Date]);
end

% Set the path to a working directory on the computer so logger data are
% transfered there and directly accessible for calculations
if TranscExtract
    WorkDir = ['~' filesep 'WorkingDirectory'];
end

%% Extracting sound events
% The samplestamp given by sound mex is not really reliable, so for each
% sound snippet, you want to find its exact location in the continuous
% recording files, then using TTL pulses, retrieve the time it correspond
% to in Deuteron, if requested.

% Checking what we have in terms of vocalization localization/extraction
ExpStartTime = DataFile(13:16);
VocExt_dir = dir(fullfile(AudioDataPath,sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)));

% Then run the logger extraction, allignment, and vocalization extraction
if TranscExtract
    fprintf(1,'*** Extract Logger data if not already done ***\n');
    % Find the ID of the recorded bats
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:P200','basic');
    RowData = find((cell2mat(RecTableData(2:end,1))== str2double(Date))) +1;
    DataInfo = RecTableData(RowData,:);
    Header = RecTableData(1,:);
    BatIDCol = find(contains(Header, 'Bat'));
    
    % extract logger data if not already done
    All_loggers_dir = dir(fullfile(Logger_dir, '*ogger*'));
    DirFlags = [All_loggers_dir.isdir];
    % Extract only those that are directories.
    All_loggers_dir = All_loggers_dir(DirFlags);
    TransceiverReset = struct(); % These are possible parameters for dealing with change of transceiver or sudden transceiver clock change. Set to empty before the first extraction
    for ll=1:length(All_loggers_dir)
        Logger_i = fullfile(Logger_dir,All_loggers_dir(ll).name);
        Ind = strfind(All_loggers_dir(ll).name, 'r');
        Logger_num = str2double(All_loggers_dir(ll).name((Ind+1):end));
        NLogCol = find(contains(Header, 'NL'));% Columns of the neural loggers
        ALogCol = find(contains(Header, 'AL'));% Columns of the audio loggers
        LogCol = NLogCol(find(cell2mat(DataInfo(NLogCol))==Logger_num));
        if isempty(LogCol) % This is an audiologger and not a neural logger
            LogCol = ALogCol(find(cell2mat(DataInfo(ALogCol))==Logger_num));
        end
        BatID = DataInfo{BatIDCol(find(BatIDCol<LogCol,1,'last'))};
        ParamFiles = dir(fullfile(Logger_i,'extracted_data','*extract_logger_data_parameters*mat'));
        if isempty(ParamFiles) || ForceExtract
            fprintf(1,'-> Extracting %s\n',All_loggers_dir(ll).name);
            
            % Bring data back on the computer
            Logger_local = fullfile(WorkDir, All_loggers_dir(ll).name);
            fprintf(1,'Transferring data from the server %s\n on the local computer %s\n', Logger_i, Logger_local);
            mkdir(Logger_local)
            [s,m,e]=copyfile(Logger_i, Logger_local, 'f');
            if ~s
                m %#ok<NOPRT>
                e %#ok<NOPRT>
                error('File transfer did not occur correctly for %s\n', Logger_i);
            end
            
            % run extraction
            if Logger_num==16 && str2double(Date)<190501
                % extract_logger_data(Logger_local, 'BatID', num2str(BatID), 'ActiveChannels', [0 1 2 3 4 5 6 7 8 9 10 12 13 14 15], 'AutoSpikeThreshFactor',5,'TransceiverReset',TransceiverReset)
                extract_logger_data(Logger_local, 'BatID', num2str(BatID), 'ActiveChannels', [0 1 2 3 4 5 6 7 8 9 10 12 13 14 15],'TransceiverReset',TransceiverReset,'AutoSpikeThreshFactor',4)
            else
                %extract_logger_data(Logger_local, 'BatID', num2str(BatID),'TransceiverReset',TransceiverReset)
                extract_logger_data(Logger_local, 'BatID', num2str(BatID),'TransceiverReset',TransceiverReset,'AutoSpikeThreshFactor',4)
            end
            
            % Keeps value of eventual clock reset
            Filename=fullfile(Logger_local, 'extracted_data', sprintf('%s_20%s_EVENTS.mat', num2str(BatID),Date));
            NewTR = load(Filename, 'TransceiverReset');
            if ~isempty(fieldnames(NewTR.TransceiverReset))% this will be used in the next loop!
                TransceiverReset = NewTR.TransceiverReset;
            end
            
            % Bring back data on the server
            fprintf(1,'Transferring data from the local computer %s\n back on the server %s\n', Logger_i, Logger_local);
            Remote_dir = fullfile(Logger_i, 'extracted_data');
            mkdir(Remote_dir)
            [s,m,e]=copyfile(fullfile(Logger_local, 'extracted_data'), Remote_dir, 'f');
            if ~s
                TicTransfer = tic;
                while toc(TicTransfer)<30*60
                    [s,m,e]=copyfile(fullfile(Logger_local, 'extracted_data'), Remote_dir, 'f');
                    if s
                        return
                    end
                end
                if ~s
                    s %#ok<NOPRT>
                    m %#ok<NOPRT>
                    e %#ok<NOPRT>
                    error('File transfer did not occur correctly for %s\n Although we tried for 30min\n', Remote_dir);
                else
                    fprintf('Extracted data transfered back on server in:\n%s\n',  Remote_dir);
                end
            else
                fprintf('Extracted data transfered back on server in:\n%s\n',  Remote_dir);
            end
            if s  %erase local data
                [sdel,mdel,edel]=rmdir(WorkDir, 's');
                if ~sdel
                    TicErase = tic;
                    while toc(TicErase)<30*60
                        [sdel,mdel,edel]=rmdir(WorkDir, 's');
                        if sdel
                            return
                        end
                    end
                end
                if ~sdel
                    sdel %#ok<NOPRT>
                    mdel %#ok<NOPRT>
                    edel %#ok<NOPRT>
                    error('File erase did not occur correctly for %s\n Although we tried for 30min\n', WorkDir);
                end
            end
            
        else
            fprintf(1,'-> Already done for %s\n',All_loggers_dir(ll).name);
        end
    end
    
    % Get the serial numbers of all audiologgers
    % bats wear
    ALThroatCol = find(contains(Header, 'AL-throat'));
    SerialNumberAL = nan(length(ALThroatCol),1);
    for dd=1:length(ALThroatCol)
        SerialNumberAL(dd) = DataInfo{ALThroatCol(dd)};
    end
    
    % Alligning TTL pulses between soundmexpro and Deuteron
    % for the RecOnly session
    
    TTL_dir = dir(fullfile(AudioDataPath,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
    if isempty(TTL_dir) || ForceAllign
        fprintf(1,'\n*** Alligning TTL pulses for the RecOnly session ***\n');
        align_soundmexAudio_2_logger(AudioDataPath, Logger_dir, ExpStartTime,'TTL_pulse_generator','Avisoft','Method','risefall', 'Session_strings', {'all voc reward stop', 'rec only stop'}, 'Logger_list', [SerialNumberAL; SerialNumberNL]);
%         align_soundmexAudio_2_logger(AudioDataPath, Logger_dir, ExpStartTime,'TTL_pulse_generator','Avisoft','Method','risefall', 'Session_strings', {'rec only start', 'rec only stop'}, 'Logger_list', [SerialNumberAL; SerialNumberNL]);
    else
        fprintf(1,'\n*** ALREADY DONE: Alligning TTL pulses for the operant session ***\n');
    end
    if isempty(VocExt_dir) || ForceVocExt1
        fprintf(1,'\n*** Localizing and extracting vocalizations that were manually extracted ***\n');
        voc_localize(VocManExtDir, AudioDataPath,Date, ExpStartTime)
    else
        fprintf(1,'\n*** ALREADY DONE: Localizing and extracting vocalizations that were manually extracted ***\n');
    end
    
    %% Identify the same vocalizations on the piezos and save sound extracts, onset and offset times
    LogVoc_dir = dir(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)));
    if isempty(LogVoc_dir) || ForceVocExt1 || ForceVocExt2
        fprintf('\n*** Localizing vocalizations on piezo recordings ***\n')
        get_logger_data_voc(AudioDataPath, Logger_dir,Date, ExpStartTime, 'SerialNumber',SerialNumberAL);
    else
        fprintf('\n*** ALREADY DONE: Localizing vocalizations on piezo recordings ***\n')
        
    end
    
    %% Identify who is calling
    fprintf('\n*** Identify who is calling ***\n')
    WhoCall_dir = dir(fullfile(Logger_dir, sprintf('*%s_%s*whocalls*', Date, ExpStartTime)));
    if isempty(WhoCall_dir) || ForceVocExt1 || ForceWhoID || ForceVocExt2
        who_calls(Logger_dir,Date, ExpStartTime,200,1);
    else
        fprintf('\n*** ALREADY DONE: Identify who is calling ***\n')
    end
    
    
elseif isempty(VocExt_dir) || ForceVocExt1
    fprintf(1,'*** Localizing and extracting vocalizations that triggered the sound detection ***\n');
    fprintf(1,'NOTE: no transceiver time extraction\n')
    voc_localize(VocManExtDir,AudioDataPath, Date, ExpStartTime,'TransceiverTime',0)
else
    fprintf(1,'*** ALREADY DONE: Localizing and extracting vocalizations that triggered the sound detection ***\n');
end





end
