function cut_neuralData_voc(Loggers_dir, Date, ExpStartTime, Flags)
load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'Voc_transc_time_refined');
NeuroBuffer = 100; % NeuroBuffer ms will be added before the onset and after the onset of the behavioral event when extracting neural data and spikes times will be alligned to behavioral event onset
MaxEventDur = NaN; % Set to NaN: The neural data is extracted for the whole duration of each event
%% Identify Neural loggers and extract the neural data that correspond to the vocalizations
% Get the number of loggers
Logger_dirs = dir(fullfile(Loggers_dir, 'logger*'));
NLogger = length(Logger_dirs);
% Identify the type of logger and extract neural data
LoggerType = cell(NLogger,1);
for ll=1:NLogger
    LData_folder = fullfile(Logger_dirs(ll).folder, Logger_dirs(ll).name,'extracted_data');
    LDir = dir(fullfile(LData_folder, '*CSC*.mat'));
    LData = load(fullfile(LDir(1).folder, LDir(1).name), 'logger_type', 'logger_serial_number');
    LoggerType{ll}  = LData.logger_type;
    if strcmp(LoggerType{ll}, 'Mous') || strcmp(LoggerType{ll}, 'Rat')
        [~,Neuro_LFP.(sprintf('Logger%s', LData.logger_serial_number)), Neuro_spikesT.(sprintf('Logger%s', LData.logger_serial_number)),Neuro_spikes.(sprintf('Logger%s', LData.logger_serial_number))] = extract_timeslot_LFP_spikes(LData_folder, Voc_transc_time_refined, NeuroBuffer,MaxEventDur, Flags);
    end
end
if Flags(2)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'Neuro_LFP', '-append');
end
if Flags(3)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)),'Neuro_spikesT', '-append');
end
if Flags(4)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)),'Neuro_spikes', '-append');
end