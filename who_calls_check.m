function [IndVocStartRaw_merge_local, IndVocStopRaw_merge_local, IndVocStartPiezo_merge_local, IndVocStopPiezo_merge_local] = who_calls_check(Raw_dir, Loggers_dir, Date, ExpStartTime, MergeThresh, Manual, varargin)

% optional parameter: Factor_RMS_Mic, Factor by which the RMS of the
% band-pass filtered baseline signal is multiplied to obtained the
% threshold of vocalization detection on Microphone

pnames = {'Working_dir','SaveFileType'};
dflts  = {Loggers_dir,'pdf'};
[Working_dir,SaveFileType] = internal.stats.parseArgs(pnames,dflts,varargin{:});
if ~exist(Working_dir,'dir')
    mkdir(Working_dir)
end
DataFile = fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime));
if ~isfile(DataFile)
    warning('Vocalization data were not extracted by get_logger_data_voc.m')
else
    load(DataFile, 'Piezo_wave', 'Piezo_FS',  'Raw_wave','FS', 'RatioRMS', 'DiffRMS','BandPassFilter', 'AudioLogs', 'RMSHigh', 'RMSLow','VocFilename');
    load(fullfile(Raw_dir, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)), 'MeanStdAmpRawExtract','Voc_filename')
    %% Identify sound elements in each vocalization extract and decide of the vocalizer
    % Output corresponds to onset and offset indices of each vocal element in
    % sound section for each vocalizing logger. A cell array of length the
    % number of sound sections identified by voc_localize or
    % voc_localize_operant. For IndVocStartRaw_merged each cell contains the onset/offset samples of each vocal element in the sound section.
    % For IndVocStartPiezo_merged each cell is a cell array of size the number of
    % audio-loggers where each cell contains the onset or offset sample of
    % individual vocal element emitted by the respective audio-loggers. Note
    % that for both these variables, vocal elements that are within MergeThresh
    % ms were merged in a single vocal element.
    %
    % Run a time window of duration 2 ms to identify who is vocalizing
    % create 2 signals for the logger, one high pass filtered above 5kHz
    % and the other one low pass filtered at 5kHz. The vocalizer will have
    % the highest energy of all vocalizers and will have more energy in the
    % lower compare to higher filtered signal
    
    
    
    if nargin<6
        Manual=0;
    end
    if nargin<5
        MergeThresh = 500; % any 2 detected call spaced by less than MergeThresh ms are merged
    end
    
    % parameters
%     Consecutive_binsMic = 10; % Number of consecutive bins of the envelope difference between highpass and low pass logger signal that has to be higher than threshold to be considered as a vocalization
    Consecutive_binsPiezo = 15; % Number of consecutive bins of the envelope difference between highpass and low pass logger signal that has to be higher than threshold to be considered as a vocalization
    Factor_RMS_low = 5.*ones(length(AudioLogs),1); % Factor by which the RMS of the low-pass filtered baseline signal is multiplied to obtained the threshold of vocalization detection on piezos
    % Factor_AmpRatio = 1.5; % Factor by which the ratio of amplitude between low and high  pass filtered baseline signals is multiplied to obtain the threshold on calling vs hearing (when the bats call there is more energy in the lower frequency band than higher frequency band of the piezo) % used to be 3
    Factor_AmpDiff = 50; % Factor by which the ratio of amplitude between low and high  pass filtered baseline signals is multiplied to obtain the threshold on calling vs hearing (when the bats call there is more energy in the lower frequency band than higher frequency band of the piezo) % used to be 3
    DB_noise = 60; % Noise threshold for the spectrogram colormap
    FHigh_spec = 90000; % Max frequency (Hz) for the raw data spectrogram
    FHigh_spec_Logger = 10000; % Max frequency (Hz) for the raw data spectrogram
    BandPassFilter = [1000 5000 9900]; % Frequency bands chosen for digital signal processing
    Fhigh_power =50; % Frequency upper bound for calculating the envelope (time running RMS)
    Fs_env = 1000; % Sample frequency of the enveloppe
    % Buffer=100; % Maximum lag calculated for the cross correlation in ms, not used in that code for identification puposes but still calculated as of now
    
    % Initialize variables
    Nvoc = length(Raw_wave);
    Amp_env_LowPassLogVoc = cell(1,Nvoc);
    Amp_env_HighPassLogVoc = cell(1,Nvoc);
    Amp_env_Mic = cell(1,Nvoc);
    LowPassLogVoc = cell(1,Nvoc);
    % LowPassMicVoc = cell(1,Nvoc);
    Fns_AL = fieldnames(Piezo_wave);
     
    % design filters of raw ambient recording, bandpass and low pass which was
    % used for the cross correlation
    [z,p,k] = butter(6,[BandPassFilter(1) 90000]/(FS/2),'bandpass');
    sos_raw_band = zp2sos(z,p,k);
    % [z,p,k] = butter(6,BandPassFilter(1:2)/(FS/2),'bandpass');
    % sos_raw_low = zp2sos(z,p,k);
    FileIn = input('Get data from the server (1) or from local working directory (2)?');
    if FileIn==1
        PreviousFile = fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, MergeThresh));
    elseif FileIn==2
        PreviousFile = fullfile(Working_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, MergeThresh));
    else
        keyboard
    end
    IndVocStartRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc)
    % a cell array of the size the number of loggers + 1 in case only one bat without a logger
    % or +2 incase no identification possible but you want to keep onset/offset of each voc and
    % for each logger the index onset of when the animal start vocalizing in the raw recording before merge
    IndVocStartPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording before merge
    IndVocStopRaw = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording before merge
    IndVocStopPiezo = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording before merge
    IndVocStart_all = cell(1,Nvoc);
    IndVocStop_all = cell(1,Nvoc);
    IndVocStartRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the raw recording
    IndVocStopRaw_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the raw recording
    IndVocStartPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index onset of when the animal start vocalizing in the piezo recording
    IndVocStopPiezo_merged = cell(1,Nvoc);% Contains for each sequence of vocalizations (Nvoc) a cell array of the size the number of loggers and for each logger the index offset of when the animal stop vocalizingin the piezo recording
    RMSRatio_all = cell(1,Nvoc);
    RMSDiff_all = cell(1,Nvoc);
    MicError = [0 0];% first element = number of corrections. second = number of detection (question)
    MicErrorType = [0 0];% first element false negative (detected as noise or already detected when it is a new call), second element false positive (vice versa)
    PiezoError = [0 0];% first element = number of corrections. second = number of detection (question)
    PiezoErrorType = [0 0]; % first element false negative (detected as noise or hearing when it is a new call), second element false positive (vice versa)
    load(PreviousFile, 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','vv','MicError','PiezoError','MicErrorType','PiezoErrorType','IndVocStartRaw','IndVocStartPiezo','IndVocStopRaw','IndVocStopPiezo');
    SaveRaw_wave=0;
    if FileIn==1 % We are in check mode, else we are restarting the check mode...
        vv=1;
    end
    %% Loop through vocalizations sequences and calculate amplitude envelopes
    for vv=vv:Nvoc
        fprintf(1,'Voc sequence %d/%d\n',vv,Nvoc);
        %% Check Raw_wave
        [Raw_wave_temp, FS] = audioread(Voc_filename{vv});
        if sum(Raw_wave{vv} - Raw_wave_temp)~=0
            warning('the sound extract is not the one expected\n')
            SaveRaw_wave=1;
            keyboard
        end
        %% Plot figure 1 with already calculated data
        % First calculate the time varying RMS of the ambient microphone
        Amp_env_LowPassLogVoc{vv} = cell(length(AudioLogs),1);
        Amp_env_HighPassLogVoc{vv} = cell(length(AudioLogs),1);
        LowPassLogVoc{vv} = cell(length(AudioLogs),1);
        
        % bandpass filter the ambient mic recording
        Filt_RawVoc = filtfilt(sos_raw_band,1,Raw_wave{vv});
        Amp_env_Mic{vv} = running_rms(Filt_RawVoc, FS, Fhigh_power, Fs_env);
        
        % Plot the spectrogram of the ambient microphone
        F1=figure(1);
        clf(F1)
        ColorCode = [get(groot,'DefaultAxesColorOrder');1 1 1; 0 1 1; 1 1 0];
        subplot(length(AudioLogs)+2,1,1)
        [Raw_Spec.to, Raw_Spec.fo, Raw_Spec.logB] = spec_only_bats(Filt_RawVoc, FS, DB_noise, FHigh_spec);
        hold on
        yyaxis right
        plot((1:length(Amp_env_Mic{vv}))/Fs_env*1000, Amp_env_Mic{vv}, 'r-', 'LineWidth',2)
        ylabel(sprintf('Amp\nMic'))
        title(sprintf('Voc %d/%d',vv,Nvoc))
        xlabel(' ')
        set(gca, 'XTick',[],'XTickLabel',{})
        
        if Manual
            pause(0.1)
            PlayerMic= audioplayer((Raw_wave{vv} - mean(Raw_wave{vv}))/std(Raw_wave{vv}), FS); %#ok<TNMLP>
            play(PlayerMic)
            pause(length(Raw_wave{vv})/FS +1)
        end
        %% Loop through the loggers, calculate envelopes and plot enveloppes and spectrograms
        for ll=1:length(AudioLogs)
            if isnan(Piezo_FS.(Fns_AL{ll})(vv)) || isempty(Piezo_wave.(Fns_AL{ll}){vv})
                fprintf(1, 'NO DATA for Vocalization %d from %s\n', vv, Fns_AL{ll})
            else
                % design the filters
                [z,p,k] = butter(6,BandPassFilter(1:2)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
                sos_low = zp2sos(z,p,k);
                [z,p,k] = butter(6,BandPassFilter(2:3)/(Piezo_FS.(Fns_AL{ll})(vv)/2),'bandpass');
                sos_high = zp2sos(z,p,k);
                % filter the loggers' signals
                if sum(isnan(Piezo_wave.(Fns_AL{ll}){vv}))~=length(Piezo_wave.(Fns_AL{ll}){vv})
                    LowPassLogVoc{vv}{ll} = (filtfilt(sos_low,1,Piezo_wave.(Fns_AL{ll}){vv})); % low-pass filter the voltage trace
                    HighPassLogVoc = (filtfilt(sos_high,1,Piezo_wave.(Fns_AL{ll}){vv})); % high-pass filter the voltage trace
                    Amp_env_LowPassLogVoc{vv}{ll}=running_rms(LowPassLogVoc{vv}{ll}, Piezo_FS.(Fns_AL{ll})(vv), Fhigh_power, Fs_env);
                    Amp_env_HighPassLogVoc{vv}{ll}=running_rms(HighPassLogVoc, Piezo_FS.(Fns_AL{ll})(vv), Fhigh_power, Fs_env);
                    
                    % Plot the low pass filtered signal of each logger
                    figure(1)
                    subplot(length(AudioLogs)+2,1,ll+1)
                    [~] = spec_only_bats(LowPassLogVoc{vv}{ll}, Piezo_FS.(Fns_AL{ll})(vv), DB_noise, FHigh_spec_Logger);
                    if ll<length(AudioLogs)
                        xlabel(' '); % supress the x label output
                        set(gca,'XTick',[],'XTickLabel',{});
                    end
                    hold on
                    yyaxis right
                    plot((1:length(Amp_env_LowPassLogVoc{vv}{ll}))/Fs_env*1000, Amp_env_LowPassLogVoc{vv}{ll}, 'b-','LineWidth', 2);
                    hold on
                    plot((1:length(Amp_env_HighPassLogVoc{vv}{ll}))/Fs_env*1000, Amp_env_HighPassLogVoc{vv}{ll}, 'r-','LineWidth',2);
                    ylabel(sprintf('Amp\n%s',Fns_AL{ll}([1:3 7:end])))
                    hold off
                    if Manual
                        pause(0.1)
                        Player= audioplayer((Piezo_wave.(Fns_AL{ll}){vv}-mean(Piezo_wave.(Fns_AL{ll}){vv}))/std(Piezo_wave.(Fns_AL{ll}){vv}), Piezo_FS.(Fns_AL{ll})(vv)); %#ok<TNMLP>
                        play(Player)
                        pause(length(Raw_wave{vv})/FS +1)
                    end
                else
                    Amp_env_LowPassLogVoc{vv}{ll}=resample(nan(1,length(Piezo_wave.(Fns_AL{ll}){vv})), Fs_env, round(Piezo_FS.(Fns_AL{ll})(vv)));
                    Amp_env_HighPassLogVoc{vv}{ll}=resample(nan(1,length(Piezo_wave.(Fns_AL{ll}){vv})), Fs_env, round(Piezo_FS.(Fns_AL{ll})(vv)));
                end
            end
        end
        if sum(cellfun('isempty',(Amp_env_LowPassLogVoc{vv}))) == length(AudioLogs) % No logger datat
            
            fprintf(1,'CANNOT DETERMINE OWNERSHIP NO DETECTION OF ONSET/OFFSET\nSet Force_Save_onoffsets_mic to 1 to force saving data\n')
            
        else % There is some data on the logger, plot the onset and offset of extracted calls
            
            for ll=1:length(AudioLogs)
                if isempty(IndVocStartRaw_merged{vv}{ll})
                    fprintf('No vocalization detected on %s\n',Fns_AL{ll});
                else% Some vocalizations were detected
                    NV = length(IndVocStartRaw_merged{vv}{ll});
                    % Now plot the onset/offset of each extract on the
                    % spectrograms of the loggers
                    figure(1)
                    for ii=1:NV
                        subplot(length(AudioLogs)+2,1,ll+1)
                        hold on
                        yyaxis right
                        YLim = get(gca, 'YLim');
                        plot([IndVocStartRaw_merged{vv}{ll}(ii) IndVocStopRaw_merged{vv}{ll}(ii)]/FS*1000, [YLim(2)*4/5 YLim(2)*4/5], 'k--', 'LineWidth',2)
                        subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
                        hold on
                        plot([IndVocStartRaw_merged{vv}{ll}(ii) IndVocStopRaw_merged{vv}{ll}(ii)]/FS*1000, [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','--')
                        hold off
                    end
                   
                end
            end
            
            figure(1)
            subplot(length(AudioLogs)+2,1,1)
            XLIM = get(gca, 'XLim');
            subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
            set(gca, 'YTick', 1:length(AudioLogs))
            set(gca, 'YTickLabel', [Fns_AL; 'Mic'])
            set(gca, 'YLim', [0 length(AudioLogs)+1])
            set(gca, 'XLim', XLIM);
            ylabel('AL ID')
            xlabel('Time (ms)')
           
        end
        ManCall = input('Do the data look fine? yes () No (anything)');
        
        if ~isempty(ManCall) % Need to re-calculate data!
            % Treat the case where some approximations of the calculations led to
            % slight diffreent sizes of running RMS
            Short = min(cellfun('length',Amp_env_LowPassLogVoc{vv}));
            Long = max(cellfun('length',Amp_env_LowPassLogVoc{vv}));
            if (Long-Short)>1
                error('The length of vectors of running RMS are too different than expected, please check!\n')
            elseif Long~=Short
                Amp_env_LowPassLogVoc{vv} = cellfun(@(X) X(1:Short), Amp_env_LowPassLogVoc{vv}, 'UniformOutput',false);
                Amp_env_HighPassLogVoc{vv} = cellfun(@(X) X(1:Short), Amp_env_HighPassLogVoc{vv}, 'UniformOutput',false);
            end
            Amp_env_LowPassLogVoc_MAT = cell2mat(Amp_env_LowPassLogVoc{vv});
            Amp_env_HighPassLogVoc_MAT = cell2mat(Amp_env_HighPassLogVoc{vv});
            RatioAmp = (Amp_env_LowPassLogVoc_MAT +1)./(Amp_env_HighPassLogVoc_MAT+1);
            DiffAmp = Amp_env_LowPassLogVoc_MAT-Amp_env_HighPassLogVoc_MAT;
            
            % Plot the ratio of time varying RMS, the difference in time varying
            % RMS between the high and low frequency bands and the absolute time
            % varying RMS of the low frequency band
            F2=figure(2);
            clf(F2)
            subplot(3,1,3)
            for ll=1:length(AudioLogs)
                plot(RatioAmp(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
                hold on
            end
            ylabel('Frequency bands Ratio')
            legend(Fns_AL{:})
            subplot(3,1,2)
            for ll=1:length(AudioLogs)
                plot(DiffAmp(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
                hold on
            end
            title('Calling vs Hearing')
            ylabel('Frequency bands Diff')
            hold on
            for ll=1:length(AudioLogs)
                plot([0 size(Amp_env_LowPassLogVoc_MAT,2)], Factor_AmpDiff * DiffRMS.(Fns_AL{ll})(1)*ones(2,1), 'Color',ColorCode(ll,:),'LineStyle','--');
                hold on
            end
            legend({Fns_AL{:} 'calling detection threshold' 'calling detection threshold'})
            subplot(3,1,1)
            for ll=1:length(AudioLogs)
                plot(Amp_env_LowPassLogVoc_MAT(ll,:), 'LineWidth',2, 'Color', ColorCode(ll,:))
                hold on
            end
            hold on
            ylabel(sprintf('Running RMS %dHz-%dHz', BandPassFilter(1:2)))
            title('Detection of vocalizations on each logger')
            for ll=1:length(AudioLogs)
                plot([0 size(Amp_env_LowPassLogVoc_MAT,2)], Factor_RMS_low(ll) * RMSLow.(Fns_AL{ll})(1)*ones(2,1), 'Color',ColorCode(ll,:),'LineStyle','--');
                hold on
            end
            legend({Fns_AL{:} 'Microphone' 'voc detection threshold' 'voc detection threshold'})
            %% Find out which calls are emitted by each individual
            Vocp = nan(size(DiffAmp) + [1 0]); % This is the output of the decision criterion to determine at each time point if a bat is vocalizing or not. each row=a logger, each column = a time point
            RowSize = length(AudioLogs);
            IndVocStartRaw{vv} = cell(RowSize,1);
            IndVocStartPiezo{vv} = cell(RowSize,1);
            IndVocStopRaw{vv} = cell(RowSize,1);
            IndVocStopPiezo{vv} = cell(RowSize,1);
            IndVocStart = cell(RowSize,1);
            IndVocStop = cell(RowSize,1);
            IndVocStartRaw_merge_local = cell(RowSize,1);
            IndVocStopRaw_merge_local = cell(RowSize,1);
            IndVocStartPiezo_merge_local = cell(RowSize,1);
            IndVocStopPiezo_merge_local = cell(RowSize,1);
            %     Xcor = cell(length(AudioLogs),1);
            RMSRatio = cell(RowSize,1);
            RMSDiff = cell(RowSize,1);
            
            for ll=1:RowSize
                Vocp(ll,:) = Amp_env_LowPassLogVoc_MAT(ll,:)>(Factor_RMS_low(ll) * RMSLow.(Fns_AL{ll})(1)); % Time points above amplitude threshold on the low-passed logger signal
                IndVocStart{ll} = strfind(Vocp(ll,:), ones(1,Consecutive_binsPiezo)); %find the first indices of every sequences of length "Consecutive_bins" higher than RMS threshold
                if isempty(IndVocStart{ll})
                    fprintf('No vocalization detected on %s\n',Fns_AL{ll});
                else% Some vocalizations were detected
                    IndVocStart_diffind = find(diff(IndVocStart{ll})>1);
                    IndVocStart{ll} = [IndVocStart{ll}(1) IndVocStart{ll}(IndVocStart_diffind +1)]; % these two lines get rid of overlapping sequences that werer detected several times
                    NV = length(IndVocStart{ll}); % This is the number of detected potential vocalization
                    IndVocStop{ll} = nan(1,NV);
                    %             Xcor{ll} = nan(NV,1);
                    RMSRatio{ll} = nan(NV,1);
                    RMSDiff{ll} = nan(NV,1);
                    Call1Hear0_temp = nan(NV,1);
                    Call1Hear0_man = nan(NV,1);
                    
                    for ii=1:NV
                        IVStop = find(Vocp(ll,IndVocStart{ll}(ii):end)==0, 1, 'first');
                        if ~isempty(IVStop)
                            IndVocStop{ll}(ii) = IndVocStart{ll}(ii) + IVStop -1;
                        else
                            IndVocStop{ll}(ii) = length(Vocp(ll,:));
                        end
                        % Plot the localization of that sound extract on figure 3
                        % Duplicate the figure of the spectrogram for manual input purposes
                        %                             if ii==1
                        F3 = figure(3);
                        clf(F3)
                        maxB = max(max(Raw_Spec.logB));
                        minB = maxB-DB_noise;
                        imagesc(Raw_Spec.to*1000,Raw_Spec.fo,Raw_Spec.logB);          % to is in seconds
                        axis xy;
                        caxis('manual');
                        caxis([minB maxB]);
                        cmap = spec_cmap();
                        colormap(cmap);
                        v_axis = axis;
                        v_axis(3)=0;
                        v_axis(4)=FHigh_spec;
                        axis(v_axis);
                        xlabel('time (ms)'), ylabel('Frequency');
                        title(sprintf('Ambient Microphone Voc %d/%d',vv,Nvoc))
                        yyaxis right
                        %                                 ylabel('Logger ID')
                        %                                 set(gca, 'YTick', 1:ll, 'YTickLabel', Fns_AL, 'YLim', [0 (length(AudioLogs)+1)],'YDir', 'reverse')
                        %                             end
                        F3 = figure(3);
                        hold on
                        yyaxis right
                        plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [ll ll], 'k:', 'LineWidth',2)
                        hold off


                        % Calculate the Average RMS Ratio for each sound extract
                        % and decide about the vocalization ownership
                        RMSRatio{ll}(ii) = mean(RatioAmp(ll,IndVocStart{ll}(ii):IndVocStop{ll}(ii)));
                        RMSDiff{ll}(ii) = mean(DiffAmp(ll,IndVocStart{ll}(ii):IndVocStop{ll}(ii)));
                        %                     Call1Hear0_temp(ii) = RMSRatio{ll}(ii) > Factor_AmpRatio * RatioRMS.(Fns_AL{ll})(1);
                        Call1Hear0_temp(ii) = RMSDiff{ll}(ii) > Factor_AmpDiff * DiffRMS.(Fns_AL{ll})(1);

                        % update figure(3) with the decision
                        figure(3)
                        yyaxis right
                        hold on
                        if Call1Hear0_temp(ii)
                            text(IndVocStop{ll}(ii)/Fs_env*1000, ll, sprintf('%s calling',Fns_AL{ll}))
                        else
                            text(IndVocStop{ll}(ii)/Fs_env*1000, ll, sprintf('%s hearing/noise',Fns_AL{ll}))
                        end
                        hold off

                        if Call1Hear0_temp(ii)
                            fprintf('Computer guess for that sound element: %s calling\n',Fns_AL{ll});
                        else
                            fprintf('Computer guess for that sound element: %s hearing/noise\n',Fns_AL{ll});
                        end
                        Call1Hear0_man(ii) = input('Indicate your choice: calling (1) hearing/noise (0) Listen again to that logger recording (any other number)\n');
                        while Call1Hear0_man(ii)~=0 && Call1Hear0_man(ii)~=1
                            Player= audioplayer((Piezo_wave.(Fns_AL{ll}){vv}-mean(Piezo_wave.(Fns_AL{ll}){vv}))/std(Piezo_wave.(Fns_AL{ll}){vv}), Piezo_FS.(Fns_AL{ll})(vv)); %#ok<TNMLP>
                            play(Player)
                            pause(length(Raw_wave{vv})/FS +1)
                            Call1Hear0_man(ii) = input('Indicate your choice: calling (1) hearing/noise (0) listen again to that logger recording (any other number)\n');
                        end
                        
                        Agree = Call1Hear0_temp(ii)== Call1Hear0_man(ii);
                        if ~Agree
                            Call1Hear0_temp(ii) = Call1Hear0_man(ii);
                            PiezoError = PiezoError + [1 1];
                            PiezoErrorType = PiezoErrorType + [Call1Hear0_temp(ii) ~Call1Hear0_temp(ii)];
                        else
                            PiezoError = PiezoError + [0 1];
                        end


                        % update figure(3) with the decision
                        figure(3)
                        yyaxis right
                        hold on
                        if Call1Hear0_temp(ii)
                            plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
                            text(IndVocStop{ll}(ii)/Fs_env*1000, ll+0.2, sprintf('%s calling',Fns_AL{ll}),'Color','r','FontWeight','bold')
                        else
                            plot([IndVocStart{ll}(ii)/Fs_env IndVocStop{ll}(ii)/Fs_env]*1000, [ll ll], 'k-', 'LineWidth',2)
                            text(IndVocStop{ll}(ii)/Fs_env*1000, ll+0.2, sprintf('%s hearing',Fns_AL{ll}),'Color','r','FontWeight','bold')
                        end
                        hold off

                    end
                        
                    % Stash noise and Only keep vocalizations that were produced by the bat
                    % according to the high RMSRatio value
                    IndVocStart{ll} = IndVocStart{ll}(logical(Call1Hear0_temp));
                    IndVocStop{ll} = IndVocStop{ll}(logical(Call1Hear0_temp));
                    IndVocStartRaw{vv}{ll} = round(IndVocStart{ll}/Fs_env*FS);
                    IndVocStartPiezo{vv}{ll} = round(IndVocStart{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                    IndVocStopRaw{vv}{ll} = round(IndVocStop{ll}/Fs_env*FS);
                    IndVocStopPiezo{vv}{ll} = round(IndVocStop{ll}/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                    NV = sum(Call1Hear0_temp);
                        
                    if NV
                        % Now merge detected vocalizations that are within MergeThresh
                        Merge01 = (IndVocStart{ll}(2:end)-IndVocStop{ll}(1:(end-1)) < (MergeThresh/1000*Fs_env));
                        NV = NV-sum(Merge01);
                        IndVocStartRaw_merge_local{ll} = nan(1,NV);
                        IndVocStopRaw_merge_local{ll} = nan(1,NV);
                        IndVocStartPiezo_merge_local{ll} = nan(1,NV);
                        IndVocStopPiezo_merge_local{ll} = nan(1,NV);
                        CutInd = find(~Merge01);
                        IndVocStartRaw_merge_local{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*FS);
                        IndVocStartPiezo_merge_local{ll}(1) = round(IndVocStart{ll}(1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                        IndVocStopRaw_merge_local{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*FS);
                        IndVocStopPiezo_merge_local{ll}(end) = round(IndVocStop{ll}(end)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                        for cc=1:length(CutInd)
                            IndVocStopRaw_merge_local{ll}(cc) = round(IndVocStop{ll}( CutInd(cc))/Fs_env*FS);
                            IndVocStopPiezo_merge_local{ll}(cc) = round(IndVocStop{ll}(CutInd(cc))/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                            IndVocStartRaw_merge_local{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*FS);
                            IndVocStartPiezo_merge_local{ll}(cc+1) = round(IndVocStart{ll}(CutInd(cc)+1)/Fs_env*Piezo_FS.(Fns_AL{ll})(vv));
                        end

                        % Now plot the onset/offset of each extract on the
                        % spectrograms of the loggers
                        figure(1)
                        for ii=1:NV
                            subplot(length(AudioLogs)+2,1,ll+1)
                            hold on
                            yyaxis right
                            YLim = get(gca, 'YLim');
                            plot([IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000, [YLim(2)*4/5 YLim(2)*4/5], 'k-', 'LineWidth',2)
                            subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
                            hold on
                            plot([IndVocStartRaw_merge_local{ll}(ii) IndVocStopRaw_merge_local{ll}(ii)]/FS*1000, [ll ll], 'Color',ColorCode(ll,:), 'LineWidth',2, 'LineStyle','-')
                            hold off
                        end
                    end
                end
            end
            figure(1)
            subplot(length(AudioLogs)+2,1,1)
            XLIM = get(gca, 'XLim');
            subplot(length(AudioLogs)+2,1,length(AudioLogs)+2)
            set(gca, 'YTick', 1:length(AudioLogs))
            set(gca, 'YTickLabel', [Fns_AL; 'Mic'])
            set(gca, 'YLim', [0 length(AudioLogs)+1])
            set(gca, 'XLim', XLIM);
            ylabel('AL ID')
            xlabel('Time (ms)')
            [~,FileVoc]=fileparts(VocFilename{vv});
            fprintf(1,'saving figures...\n')
            if strcmp(SaveFileType,'pdf')
                print(F1,fullfile(Working_dir,sprintf('%s_whocalls_spec_%d.pdf', FileVoc, MergeThresh)),'-dpdf','-fillpage')
                if ManCall % Only save the RMS figure if there was a vocalization
                    saveas(F2,fullfile(Working_dir,sprintf('%s_whocalls_RMS_%d.pdf', FileVoc, MergeThresh)),'pdf')
                end
            elseif strcmp(SaveFileType,'fig')
                saveas(F1,fullfile(Working_dir,sprintf('%s_whocalls_spec_%d.fig', FileVoc, MergeThresh)))
                if ManCall % Only save the RMS figure if there was a vocalization
                    saveas(F2,fullfile(Working_dir,sprintf('%s_whocalls_RMS_%d.fig', FileVoc, MergeThresh)))
                end
            end
            if ManCall
                pause(1)
            end
            clf(F1)
            clf(F2)
            
            % Gather Vocalization production data:
            IndVocStart_all{vv} = IndVocStart;
            IndVocStop_all{vv} = IndVocStop;
            IndVocStartRaw_merged{vv} = IndVocStartRaw_merge_local;
            IndVocStopRaw_merged{vv} = IndVocStopRaw_merge_local;
            IndVocStartPiezo_merged{vv} = IndVocStartPiezo_merge_local;
            IndVocStopPiezo_merged{vv} = IndVocStopPiezo_merge_local;
            RMSRatio_all{vv} = RMSRatio;
            RMSDiff_all{vv} = RMSDiff;
            clf(F2)
        end
        fprintf(1,'saving data...\n')
        save(fullfile(Working_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'IndVocStartRaw', 'IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo', 'IndVocStart_all', 'IndVocStop_all','RMSRatio_all','RMSDiff_all','vv','MicError','PiezoError','MicErrorType','PiezoErrorType');
        clf(F1)
        
    end
    if ~strcmp(Working_dir,Loggers_dir)
        fprintf(1,'Transferring data back on the server\n')
        [s,m,e]=copyfile(fullfile(Working_dir,'*'), Loggers_dir, 'f');
        if ~s
            fprintf(1,'File transfer did not occur correctly\n')
            keyboard
        end
        if s  %erase local data
            [sdel,mdel,edel]=rmdir(Working_dir, 's');
        end
    end
    if SaveRaw_wave
        save(DataFile, 'Raw_wave','-append')
    end
end
end