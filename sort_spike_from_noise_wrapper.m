function sort_spike_from_noise_wrapper(Loggers_dir,RThresh,PLOT)
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
        fprintf(1,' -> Logger %s\n', LData.logger_serial_number)
        % List Spike time files (unsorted tetrode spikes)
        ST_files = dir(fullfile(LData_folder, '*Tetrode_spikes_snippets*.mat'));
        NT = length(ST_files);
        % loop through them and calculate spike correlations
        parfor tt=1:NT
            fprintf(1,'   ... Loading T%d\n', tt)
            sort_spike_from_noise(fullfile(ST_files(tt).folder, ST_files(tt).name),RThresh,PLOT)
        end
    end
end
end