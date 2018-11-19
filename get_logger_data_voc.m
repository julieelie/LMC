function get_logger_data_voc(Audio_dir, Loggers_dir, Date, ExpStartTime)

Buffer = 100; % Let's cut the audio extracts Buffer ms before and after the predicted time according to audio/transceiver allignment to better allign
BandPassFilter = [1000 5000 9000];
%Parameter for detecting who is vocalizing:
Fhigh_power = 20; %Hz
Fs_env = 1000; %Hz Sample frequency of the enveloppe
Dur_RMS = 10; % duration of the sample in min for the calculation of average running RMS
%% Load the localization info of vocalization extracts
VocExt=load(fullfile(Audio_dir, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)));
% Get the number of vocalizations
Nvoc = length(VocExt.Voc_filename);
fprintf('****** Realligning and extracting data for %d vocalizations. ******\n', Nvoc);

%% Get loggers info and initialize output variables
% Get the number of loggers
Logger_dirs = dir(fullfile(Loggers_dir, 'logger*'));
NLogger = length(Logger_dirs);
% Identify the type of logger and initialize output variables
LoggerType = cell(NLogger,1);
for ll=1:NLogger
    LDir = dir(fullfile(Logger_dirs(ll).folder, Logger_dirs(ll).name, 'extracted_data', '*CSC*.mat'));
    LData = load(fullfile(LDir(1).folder, LDir(1).name), 'logger_type', 'logger_serial_number');
    LoggerType{ll}  = LData.logger_type;
    if strcmp(LoggerType{ll}, 'Audi')
        Piezo_wave.(sprintf('Logger%s', LData.logger_serial_number)) = cell(Nvoc,1);
        Piezo_FS.(sprintf('Logger%s', LData.logger_serial_number)) = nan(Nvoc,1);
    end
end
Raw_wave = cell(Nvoc,1);

%% Extract the audio loggers data that correspond to the vocalizations
fprintf(1, '*** Extract audio loggers data for each vocalization ***\n');
AudioLogs = find(contains(LoggerType, 'Audi'));
for ll=1:length(AudioLogs)
    fprintf(1, 'AL %d/%d\n', ll, length(AudioLogs));
    % load the audio data
    LDir = dir(fullfile(Logger_dirs(AudioLogs(ll)).folder, Logger_dirs(AudioLogs(ll)).name,'extracted_data', '*CSC*.mat'));
    if length(LDir)~=1
        error('There are %d CSC files when there should be only one for this audio logger: %s', length(LDir),Logger_dirs(AudioLogs(ll)).name)
    end
    LData = load(fullfile(LDir(1).folder, LDir(1).name));
    % Get the average running rms in the lower frequency band and higher
    % frequency band for that logger in a 5 min extract in the middle of
    % the recording
    fprintf(1, 'Calculating average RMS values on a %d min sample\n',Dur_RMS);
    FS_Logger_local = nanmean(LData.Estimated_channelFS_Transceiver);
    SampleDur = round(Dur_RMS*60*FS_Logger_local);
    StartSamp = round(length(LData.AD_count_int16)/2);
    [z,p,k] = butter(6,BandPassFilter(1:2)/(round(FS_Logger_local)/2),'bandpass');
    sos_low = zp2sos(z,p,k);
    [z,p,k] = butter(6,BandPassFilter(2:3)/(round(FS_Logger_local)/2),'bandpass');
    sos_high = zp2sos(z,p,k);
    Filtered_voltage_trace = filtfilt(sos_low,1,double(LData.AD_count_int16(StartSamp + (1:SampleDur))));
    Amp_env_voltage_low=running_rms(Filtered_voltage_trace, FS_Logger_local, Fhigh_power, Fs_env);
    
    Filtered_voltage_trace = filtfilt(sos_high,1,double(LData.AD_count_int16(StartSamp + (1:SampleDur))));
    Amp_env_voltage_high=running_rms(Filtered_voltage_trace, FS_Logger_local, Fhigh_power, Fs_env);
    RMSHigh.(sprintf('Logger%s', LData.logger_serial_number))(1) = mean(Amp_env_voltage_high);
    RMSHigh.(sprintf('Logger%s', LData.logger_serial_number))(2) = std(Amp_env_voltage_high);
    RMSLow.(sprintf('Logger%s', LData.logger_serial_number))(1) = mean(Amp_env_voltage_low);
    RMSLow.(sprintf('Logger%s', LData.logger_serial_number))(2) = std(Amp_env_voltage_low);
    RatioRMS.(sprintf('Logger%s', LData.logger_serial_number))(1) = mean(Amp_env_voltage_low ./ Amp_env_voltage_high);
    RatioRMS.(sprintf('Logger%s', LData.logger_serial_number))(2) = std(Amp_env_voltage_low ./ Amp_env_voltage_high);
    DiffRMS.(sprintf('Logger%s', LData.logger_serial_number))(1) = mean(Amp_env_voltage_low - Amp_env_voltage_high);
    DiffRMS.(sprintf('Logger%s', LData.logger_serial_number))(2) = std(Amp_env_voltage_low - Amp_env_voltage_high);
    clear Amp_env_voltage_low Amp_env_voltage_high Filtered_voltage_trace
    
    % loop through vocalizations and extract the snippet of logger and
    % audio data
    for vv=1:Nvoc
        fprintf(1, 'Extracting vocalization %d/%d\n', vv, Nvoc);
        VocOnset_time = VocExt.Voc_transc_time(vv,1) - Buffer;
        VocOffset_time = VocExt.Voc_transc_time(vv,2) + Buffer;
        % find the time stamp on the logger that is closest to before
        % the snippet of sound onset
        IndTSOn = find(LData.Timestamps_of_first_samples_usec<(VocOnset_time*10^3), 1, 'Last');
        % find the time stamp on the logger that is closest to after
        % the snippet of sound offset
        IndTSOff = find(LData.Timestamps_of_first_samples_usec>(VocOffset_time*10^3), 1, 'First');
        % deduct the corresponding onset and offset samples
        IndSampOn = round(LData.Indices_of_first_and_last_samples(IndTSOn,1) + LData.Estimated_channelFS_Transceiver(IndTSOn)*(10^-6)*(VocOnset_time*10^3 - LData.Timestamps_of_first_samples_usec(IndTSOn)));
        IndSampOff = round(LData.Indices_of_first_and_last_samples(IndTSOff,1) - LData.Estimated_channelFS_Transceiver(IndTSOff)*(10^-6)*(LData.Timestamps_of_first_samples_usec(IndTSOff) - VocOffset_time*10^3));
        % extract the data snippet
%         Piezo_wave.(sprintf('Logger%s', LData.logger_serial_number)){vv} = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16))/std(LData.AD_count_int16);
        Piezo_wave.(sprintf('Logger%s', LData.logger_serial_number)){vv} = double(LData.AD_count_int16(IndSampOn:IndSampOff));
        Piezo_FS.(sprintf('Logger%s', LData.logger_serial_number))(vv) = mean(LData.Estimated_channelFS_Transceiver(IndTSOn:IndTSOff));
    end
end
clear LData
 
%% Roughly Identify who is vocalizing and re-adjust the allignement by doing a cross-correlation
fprintf(1, 'Ajusting allignment of Microphone data to audio loggers data for each vocalization\n');
Voc_transc_time_refined = nan(Nvoc,2); % New more accurate time in transceiver time of vocalization onset/offset on loggers in ms
OnsetAudiosamp = nan(Nvoc,1); % New more accurate onset sample of vocalization on audio logger extracts in Piezo_wave
OffsetAudiosamp=nan(Nvoc,1); % New more accurate offset sample of vocalization on audio logger extracts in Piezo_wave
Fns_AL = fieldnames(Piezo_wave);

for vv=1:Nvoc
    % calculate the RMS of each logger and select the one with the highest
    % to perform a cross-correlation
    RMS_audioLog = nan(length(AudioLogs),1);
    Filt_Logger_wav = cell(length(AudioLogs),1);
    F1=figure(1);
    for ll=1:length(AudioLogs)
        [z,p,k] = butter(6,BandPassFilter(1:2)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
        sos = zp2sos(z,p,k);
        Filt_Logger_wav{ll} = (filtfilt(sos,1,Piezo_wave.(Fns_AL{ll}){vv} - mean(Piezo_wave.(Fns_AL{ll}){vv}))); % band-pass filter the centered voltage trace
        RMS_audioLog(ll) = mean(Filt_Logger_wav{ll} .^2)^0.5;
        subplot(length(AudioLogs)+1,1,ll)
        plot(Piezo_wave.(Fns_AL{ll}){vv} - mean(Piezo_wave.(Fns_AL{ll}){vv}), 'k-');
        hold on
        plot(Filt_Logger_wav{ll}, 'r-')
        if ll==1
            title('Loggers signal filtering and resampling')
            legend('centered voltage trace', sprintf('BandPass %d %d Hz', BandPassFilter(1:2)))
        end
        hold off
    end
    % Only rely on the logger of the individual that vocalized the loudest
    % for reallignement
    HRMS_Ind = find(RMS_audioLog == max(RMS_audioLog));
    
    % filter the original wavfile
    [Raw_wave{vv}, FS] = audioread(VocExt.Voc_filename{vv});
    [z,p,k] = butter(6,BandPassFilter(1:2)/(FS/2),'bandpass');% a 12th order Butterworth band-pass filter; the second input argument is normalized cut-off frequency (ie. normalized to the Nyquist frequency, which is half the sampling frequency, as required by MATLAB)
    sos = zp2sos(z,p,k); % obtain the second order section (biquad) filter to use as input in filtfilt
    Filt_Raw_wav=(filtfilt(sos,1,Raw_wave{vv})); % band-pass filter the voltage traces
    subplot(length(AudioLogs)+1,1,length(AudioLogs)+1)
    plot(Raw_wave{vv}, '-k')
    hold on
    plot(Filt_Raw_wav, '-r')
    hold off
    title('Environmental Mic filtering and resampling')
    legend('Raw voltage trace', sprintf('BandPass %d %d Hz', BandPassFilter(1:2)))

    % resample the sounds so they are at the same sample frequency of 4
    % times the low pass filter value
    Resamp_Filt_Raw_wav = resample(Filt_Raw_wav, 4*BandPassFilter(2), FS);
    subplot(length(AudioLogs)+1,1,length(AudioLogs)+1)
    t2 = (0:(length(Resamp_Filt_Raw_wav)-1))*FS/(4*BandPassFilter(2));
    hold on
    plot(t2, Resamp_Filt_Raw_wav, 'g-')
    legend(sprintf('BandPass + Resampled %d Hz', 4*BandPassFilter(2)))
    hold off
    Resamp_Filt_Logger_wav = cell(length(HRMS_Ind),1);
    for ll=1:length(HRMS_Ind)
        FS_ll = round(Piezo_FS.(Fns_AL{HRMS_Ind(ll)})(vv));
        Resamp_Filt_Logger_wav{ll} = resample(Filt_Logger_wav{HRMS_Ind(ll)}, 4*BandPassFilter(2),FS_ll);
        subplot(length(AudioLogs)+1,1,HRMS_Ind(ll))
        t2 = (0:(length(Resamp_Filt_Logger_wav{ll})-1))*FS_ll/(4*BandPassFilter(2));
        hold on
        plot(t2, Resamp_Filt_Logger_wav{ll} , 'g-')
        legend(sprintf('BandPass + Resampled %d Hz', 4*BandPassFilter(2)))
        hold off
    end
    
    % Localize vocalizations on the selected loggers and set to zero the
    % signal outside of putative vocalizations.
    % calculate a running RMS of the audio signal with a 1ms bin
    % resolution. a vocalization is defined as 50 consecutive time bins
    % (50ms) above the mean running RMS here
    Consecutive_bins = 50;
    Fhigh_power = 20; %Hz
    Fs_env = 1000; %Hz
    Amp_env_voltage = cell(length(HRMS_Ind),1);
    Clean_Resamp_Filt_Logger_wav = cell(length(HRMS_Ind),1);
    for ll=1:length(HRMS_Ind)
        Amp_env_voltage{ll}=running_rms(Resamp_Filt_Logger_wav{ll}, 4*BandPassFilter(2), Fhigh_power, Fs_env);
        Vocp = Amp_env_voltage{ll}>median(Amp_env_voltage{ll});
        IndVocStart = strfind(Vocp, ones(1,Consecutive_bins));
        IndVocStart_diffind = find(diff(IndVocStart)>1);
        IndVocStart = [IndVocStart(1) IndVocStart(IndVocStart_diffind +1)];
        NV = length(IndVocStart);
        IndVocStop = nan(NV,1);
        Clean_Resamp_Filt_Logger_wav{ll} = zeros(size(Resamp_Filt_Logger_wav{ll}));
        for ii=1:NV
            IVStop = find(Vocp(IndVocStart(ii):end)==0, 1, 'first');
            if isempty(IVStop)
                IVStop = length(Vocp(IndVocStart(ii):end));
            end
            IndVocStop(ii) = IndVocStart(ii) + IVStop;
            IndVocStart(ii) = round(IndVocStart(ii)/Fs_env*4*BandPassFilter(2));
            IndVocStop(ii) = round(IndVocStop(ii)/Fs_env*4*BandPassFilter(2));
            if IndVocStop(ii)>length(Clean_Resamp_Filt_Logger_wav{ll}) % This sound element reaches the end of the recording and downsampling messed up the exact indices values
                Clean_Resamp_Filt_Logger_wav{ll}(IndVocStart(ii):end) = Resamp_Filt_Logger_wav{ll}(IndVocStart(ii):end);
            else
                Clean_Resamp_Filt_Logger_wav{ll}(IndVocStart(ii):IndVocStop(ii)) = Resamp_Filt_Logger_wav{ll}(IndVocStart(ii):IndVocStop(ii));
            end
        end
        figure(1)
        subplot(length(AudioLogs)+1,1,HRMS_Ind(ll))
        t2 = (0:(length(Resamp_Filt_Logger_wav{ll})-1))*FS_ll/(4*BandPassFilter(2));
        hold on
        plot(t2, Clean_Resamp_Filt_Logger_wav{ll} , 'c-')
        legend('Cleaned signal for cross-correlation')
        hold off
    end
    
    
    % cross correlate the raw data with the logger and detect what the lag
    % is.
    Xcor = cell(length(HRMS_Ind),1);
    Lag = cell(length(HRMS_Ind),1);
    LagDiff = nan(length(HRMS_Ind),1);
    F2 = figure(2);
    for ll=1:length(HRMS_Ind)
        [Xcor{ll},Lag{ll}] = xcorr(Resamp_Filt_Raw_wav,Clean_Resamp_Filt_Logger_wav{ll}, (Buffer*2)*4*BandPassFilter(2)); % Running a cross correlation between the raw signal and each audio logger signal with a maximum lag equal to twice the Buffer size
        [~,I] = max(abs(Xcor{ll}));
        LagDiff(ll) = Lag{ll}(I);
        subplot(length(HRMS_Ind),1,ll)
        plot(Lag{ll},Xcor{ll})
        hold on
        line([LagDiff(ll) LagDiff(ll)], [min(Xcor{ll}) max(Xcor{ll})], 'Color','r', 'LineStyle','--');
        text(LagDiff(ll) , mean(Xcor{ll}), sprintf('%d',LagDiff(ll)))
        hold off
    end
    
    % check the allignment
    F3=figure(3);
    subplot(length(HRMS_Ind)+1,1,1)
    plot(Resamp_Filt_Raw_wav, 'k-')
    title(sprintf('Mic filtered resampled data Voc %d/%d',vv,Nvoc))
    for ll=1:length(HRMS_Ind)
        subplot(length(HRMS_Ind)+1,1,ll+1)
        plot(Resamp_Filt_Logger_wav{ll}(-LagDiff(ll) + (0:(length(Resamp_Filt_Logger_wav{ll}) - (Buffer*2*10^-3)*4*BandPassFilter(2)))), 'k')
        hold on
        plot(Clean_Resamp_Filt_Logger_wav{ll}(-LagDiff(ll) + (0:(length(Clean_Resamp_Filt_Logger_wav{ll}) - (Buffer*2*10^-3)*4*BandPassFilter(2)))), 'c')
        title(sprintf('Logger filtered resampled data Voc %d/%d',vv,Nvoc))
    end
    pause(1)
    Player= audioplayer((Resamp_Filt_Raw_wav -mean(Resamp_Filt_Raw_wav))/std(Resamp_Filt_Raw_wav), 4*BandPassFilter(2)); %#ok<TNMLP>
    play(Player)
    pause()
    clf(F1)
    clf(F2)
    clf(F3)
    LagDiff = mean(LagDiff);
    
    % calculating the portion of data to erase or add in front of the
    % vocalization on the logger recordings.
    TimeDiff_audio = LagDiff/(4*BandPassFilter(2)); % this should be negative the reference being the vocalization onset -buffer of 100ms
    OnsetAudiosamp(vv) = round(-TimeDiff_audio*Piezo_FS.(Fns_AL{ll})(vv));
    OffsetAudiosamp(vv) = round(length(Piezo_wave.(Fns_AL{ll}){vv}) - (2*Buffer*(10^-3)*Piezo_FS.(Fns_AL{ll})(vv) - OnsetAudiosamp(vv)));
    Voc_transc_time_refined(vv,1) = Buffer + LagDiff/(4*BandPassFilter(2)) + VocExt.Voc_transc_time(vv,1); % This is the new estimate of VocExt.Voc_transc_time(vv,1) in ms
    Voc_transc_time_refined(vv,2) = Voc_transc_time_refined(vv,1) + diff(VocExt.Voc_transc_time(vv,:)); % This is the new estimate of VocExt.Voc_transc_time(vv,2) in ms
    
end

%% Better cut the audio loggers data that correspond to the vocalizations
fprintf(1, 'Better cut audio loggers data\n') 
for ll=1:length(AudioLogs)
fprintf(1,'-Audio %s\n',Fns_AL{ll})
        for vv=1:Nvoc
            Piezo_wave.(Fns_AL{ll}){vv} = Piezo_wave.(Fns_AL{ll}){vv}(OnsetAudiosamp(vv):OffsetAudiosamp(vv));
        end
end
% Save data
VocFilename = VocExt.Voc_filename;
save(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'Piezo_wave', 'Piezo_FS',  'Raw_wave','FS', 'RatioRMS', 'DiffRMS','BandPassFilter', 'AudioLogs', 'RMSHigh', 'RMSLow','VocFilename','Voc_transc_time_refined','LoggerType');

    


