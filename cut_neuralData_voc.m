function cut_neuralData_voc(Loggers_dir, Date, ExpStartTime, Flags, NeuroBuffer)
%% This function uses the better estimation of vocalization onset/offset in transceiver time calculated by get_logger_data_voc
% (Voc_transc_time_refined) And extract the corresponding neural data in the neural loggers
load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'Voc_transc_time_refined');
if nargin<5
    NeuroBuffer = 100; % NeuroBuffer ms will be added before the onset and after the offset of the behavioral event when extracting neural data and spikes times will be alligned to behavioral event onset
end
MaxEventDur = NaN; % Set to NaN: The neural data is extracted for the whole duration of each event
%% Identify Neural loggers and extract the neural data that correspond to the vocalizations
% Get the number of loggers
Logger_dirs = dir(fullfile(Loggers_dir, '*ogger*'));
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
        [~,Neuro_LFP.(sprintf('Logger%s', LData.logger_serial_number)), Neuro_spikesT.(sprintf('Logger%s', LData.logger_serial_number)),Neuro_spikes.(sprintf('Logger%s', LData.logger_serial_number))] = extract_timeslot_LFP_spikes(LData_folder, Voc_transc_time_refined, NeuroBuffer,MaxEventDur, Flags);
    end
end
if Flags(2)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)), 'Neuro_LFP', '-append');
end
if Flags(3)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)),'Neuro_spikesT', '-append');
end
if Flags(4)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)),'Neuro_spikes', '-append');
end