function cut_neuralData_behav(Loggers_dir, Date, ExpStartTime, Flags, NeuroBuffer, BehaviorType)
%% This function uses the estimation of behavior onset/offset in transceiver time calculated by get_logger_data_behav
% (AllActions_Time) And extract the corresponding neural data in the neural loggers
load(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData.mat', Date, ExpStartTime)), 'AllActions_Time', 'AllActions_ID','UActionText');
if nargin<5
    NeuroBuffer = 100; % NeuroBuffer ms will be added before the onset and after the offset of the behavioral event when extracting neural data and spikes times will be alligned to behavioral event onset
end
if nargin<6
    BehaviorType = {'licking' 'chewing' 'quiet' 'teeth-cleaning'};
end
MaxEventDur = NaN; % Set to NaN: The neural data is extracted for the whole duration of each event
%% Identify Neural loggers and extract the neural data that correspond to the Behaviors
% Get the number of loggers
Logger_dirs = dir(fullfile(Loggers_dir, '*ogger*'));
Logger_dirs=Logger_dirs([Logger_dirs.isdir]);
NLogger = length(Logger_dirs);
% Loop through behavior type and Logger
% Identify the type of logger and extract neural data
LoggerType = cell(NLogger,1);
for ll=1:NLogger
    LData_folder = fullfile(Logger_dirs(ll).folder, Logger_dirs(ll).name,'extracted_data');
    LDir = dir(fullfile(LData_folder, '*CSC*.mat'));
    LData = load(fullfile(LDir(1).folder, LDir(1).name), 'logger_type', 'logger_serial_number');
    LoggerType{ll}  = LData.logger_type;
    if strcmp(LoggerType{ll}, 'Mous') || strcmp(LoggerType{ll}, 'Rat')
        if Flags(4)    
            SUFiles = dir(fullfile(LData_folder, '*TT*SS*.mat'));
            if isempty(SUFiles)
               SUNTTFiles = dir(fullfile(LData_folder, '*TT*SS*.ntt'));
               if isempty(SUNTTFiles)
                   error('No spike sorted files available although these data are requested\n');
               else
                   error('ntt files of spike sorted data need to be converted to mat files on a windows machine!\n');
               end
            end
        end
        % Loop thourgh behavior type
        Neuro_Raw.(sprintf('Logger%s', LData.logger_serial_number)) = cell(length(UActionText),1);
        Neuro_LFP.(sprintf('Logger%s', LData.logger_serial_number)) = cell(length(UActionText),1);
        Neuro_spikesT.(sprintf('Logger%s', LData.logger_serial_number)) = cell(length(UActionText),1);
        Neuro_spikes.(sprintf('Logger%s', LData.logger_serial_number)) = cell(length(UActionText),1);
        for bb=1:length(BehaviorType)
            IndBehav = strcmp(UActionText, BehaviorType(bb));
            [Neuro_Raw.(sprintf('Logger%s', LData.logger_serial_number)){IndBehav}, Neuro_LFP.(sprintf('Logger%s', LData.logger_serial_number)){IndBehav}, Neuro_spikesT.(sprintf('Logger%s', LData.logger_serial_number)){IndBehav},Neuro_spikes.(sprintf('Logger%s', LData.logger_serial_number)){IndBehav}] = extract_timeslot_LFP_spikes(LData_folder, AllActions_Time{IndBehav}, NeuroBuffer,MaxEventDur, Flags);
        end
    end
end
if Flags(1)
    save(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)), 'Neuro_Raw', '-append');
end
if Flags(2)
    save(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)), 'Neuro_LFP', '-append');
end
if Flags(3)
    save(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)),'Neuro_spikesT', '-append');
end
if Flags(4)
    save(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)),'Neuro_spikes', '-append');
end

% Plot PSTH for each behavior for each neural logger 16
for ll=1:NLogger
    LData_folder = fullfile(Logger_dirs(ll).folder, Logger_dirs(ll).name,'extracted_data');
    LDir = dir(fullfile(LData_folder, '*CSC*.mat'));
    LData = load(fullfile(LDir(1).folder, LDir(1).name), 'logger_type', 'logger_serial_number');
    LoggerType{ll}  = LData.logger_type;
    if strcmp(LoggerType{ll}, 'Mous') || strcmp(LoggerType{ll}, 'Rat')
        plot_psth_one_voc(Neuro_spikesT.(sprintf('Logger%s', LData.logger_serial_number)), Raw_wave, Piezo_wave.Logger6, Piezo_FS.Logger6, FS)
    end
end