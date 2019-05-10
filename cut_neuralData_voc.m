function cut_neuralData_voc(Loggers_dir, Date, ExpStartTime, Flags, NeuroBuffer, DenoiseT, Rthreshold)
%% This function uses the better estimation of vocalization onset/offset in transceiver time (ms) calculated by get_logger_data_voc
% (Voc_transc_time_refined) And extract the corresponding neural data in
% the neural loggers as long as neural data of baseline activity in the
% seconds preceding the vocalization
load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'Voc_transc_time_refined', 'Raw_wave','Piezo_wave', 'FS', 'Piezo_FS');
if nargin<5
    NeuroBuffer = 100; % NeuroBuffer ms will be added before the onset and after the offset of the behavioral event when extracting neural data and spikes times will be alligned to behavioral event onset
end
if nargin<6
    DenoiseT = 0; % No sort of the tetrode spike from noise
end
if nargin<7
    Rthreshold = [0.92 0.94 0.96 0.98];
end
MaxEventDur = NaN; % Set to NaN: The neural data is extracted for the whole duration of each event
BaselineDur = 1000; % Duration of the baseline section that is seeked at least 1 second before sound onset
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
        LocalFlags = Flags;
        if Flags(4)    
            SUFiles = dir(fullfile(LData_folder, '*TT*SS*.mat'));
            if isempty(SUFiles)
               SUNTTFiles = dir(fullfile(LData_folder, '*TT*SS*.ntt'));
               if isempty(SUNTTFiles)
                   warning('No spike sorted files available for %s although these data are requested\nOnly Tetrode data will be extracted\n',Logger_dirs(ll).name);
                   LocalFlags(4) = 0;
               else
                   error('ntt files of spike sorted data need to be converted to mat files on a windows machine!\n');
               end
            end
        end
        [Neuro_Raw.(sprintf('Logger%s', LData.logger_serial_number)), Neuro_LFP.(sprintf('Logger%s', LData.logger_serial_number)), Neuro_spikesT.(sprintf('Logger%s', LData.logger_serial_number)),Neuro_spikes.(sprintf('Logger%s', LData.logger_serial_number)),Neuro_spikesTDeNoiseInd.(sprintf('Logger%s', LData.logger_serial_number))] = extract_timeslot_LFP_spikes(LData_folder, Voc_transc_time_refined, NeuroBuffer,MaxEventDur, LocalFlags, DenoiseT,Rthreshold);       
    end
end
if Flags(1)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)), 'Neuro_Raw', '-append');
end
if Flags(2)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)), 'Neuro_LFP', '-append');
end
if Flags(3)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)),'Neuro_spikesT','Neuro_spikesTDeNoiseInd', '-append');
end
if Flags(4)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)),'Neuro_spikes', '-append');
end

% % Plot PSTH for each vocalization for each neural logger 16
% for ll=1:NLogger
%     LData_folder = fullfile(Logger_dirs(ll).folder, Logger_dirs(ll).name,'extracted_data');
%     LDir = dir(fullfile(LData_folder, '*CSC*.mat'));
%     LData = load(fullfile(LDir(1).folder, LDir(1).name), 'logger_type', 'logger_serial_number');
%     LoggerType{ll}  = LData.logger_type;
%     if strcmp(LoggerType{ll}, 'Mous') || strcmp(LoggerType{ll}, 'Rat')
%         plot_psth_one_voc(Neuro_spikesT.(sprintf('Logger%s', LData.logger_serial_number)), Raw_wave, Piezo_wave.Logger6, Piezo_FS.Logger6, FS)
%     end
% end

%% Identify Neural loggers and extract the neural data that correspond to silence the second before the vocalizations

% Find sections of 1 seconds before 1s of each vocalization event
NVoc = size(Voc_transc_time_refined,1);
BSL_transc_time_refined = nan(size(Voc_transc_time_refined));
for vv=1:NVoc
    PotentialOnset = Voc_transc_time_refined(vv,1)-1000-BaselineDur;
    if vv==1 % Nothing to fear before
        BSL_transc_time_refined(vv,:) = PotentialOnset + [0 BaselineDur];
    else
        if PotentialOnset>Voc_transc_time_refined(vv-1,2) % all good to go
            BSL_transc_time_refined(vv,:) = PotentialOnset + [0 BaselineDur];
        else %take at least one second after previous vocalization offset and up to one second before current vocalization onset
            BSL_transc_time_refined(vv,1) = Voc_transc_time_refined(vv-1,2) + 1000;
            BSL_transc_time_refined(vv,2) = Voc_transc_time_refined(vv,1) - 1000;
            if diff(BSL_transc_time_refined(vv,:))<0 % No baseline available here!
                BSL_transc_time_refined(vv,:) = nan(1,2);
            end
                
        end
    end
end

% Identify the type of logger and extract neural data
IndLog = find(contains(LoggerType, 'Mous'));
for ll=1:length(IndLog)
    LData_folder = fullfile(Logger_dirs(IndLog(ll)).folder, Logger_dirs(IndLog(ll)).name,'extracted_data');
    LDir = dir(fullfile(LData_folder, '*CSC*.mat'));
    LData = load(fullfile(LDir(1).folder, LDir(1).name), 'logger_type', 'logger_serial_number');
    LocalFlags = Flags;
    if Flags(4)
        SUFiles = dir(fullfile(LData_folder, '*TT*SS*.mat'));
        if isempty(SUFiles)
            SUNTTFiles = dir(fullfile(LData_folder, '*TT*SS*.ntt'));
            if isempty(SUNTTFiles)
                warning('No spike sorted files available for %s although these data are requested\nOnly Tetrode data will be extracted\n',Logger_dirs(ll).name);
                LocalFlags(4) = 0;
            else
                error('ntt files of spike sorted data need to be converted to mat files on a windows machine!\n');
            end
        end
    end
    [Neuro_Raw_Baseline.(sprintf('Logger%s', LData.logger_serial_number)), Neuro_LFP_Baseline.(sprintf('Logger%s', LData.logger_serial_number)), Neuro_spikesT_Baseline.(sprintf('Logger%s', LData.logger_serial_number)),Neuro_spikes_Baseline.(sprintf('Logger%s', LData.logger_serial_number))] = extract_timeslot_LFP_spikes(LData_folder, BSL_transc_time_refined, 0,MaxEventDur, LocalFlags);       
end
if Flags(1)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)), 'Neuro_Raw_Baseline', '-append');
end
if Flags(2)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)), 'Neuro_LFP_Baseline', '-append');
end
if Flags(3)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)),'Neuro_spikesT_Baseline', '-append');
end
if Flags(4)
    save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)),'Neuro_spikes_Baseline', '-append');
end
save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, NeuroBuffer)),'BSL_transc_time_refined', '-append');
end