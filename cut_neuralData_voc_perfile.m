function cut_neuralData_voc_perfile(InputDataFile, OutputPath, NeuroBuffer, Flags, DenoiseT, Rthreshold)
%% This function uses the better estimation of vocalization boouts onset/offset in transceiver time (ms) calculated by get_logger_data_voc
% (Voc_transc_time_refined) And extract the corresponding neural data in
% the input datafile. A Buffer in ms is used to silghtly enlarged the
% window given by Voc_transc_time_refined. It is set to be the buffer used
% to merge calls in bouts when these are parsed out and individualized by who_calls.
% This is reasonnable here to assume that there is no vocalizations if
% NeuroBuffer is set to be 200ms, but is dangerous for the first and last
% vocalizations of the bouts as NeuroBuffer becomes larger.
% This functions also extracts neural data of baseline activity in the
% seconds preceding the vocalization, this baseline calculation assumes
% that there is no vocalization or behavioral event between the behavioral
% events indicated by Voc_trans_time_refined.
% 

% INPUT:
%       InputDataFile: path to a tetrode mat file, a single unit mat file or the
%       raw data (*CSC*.mat file)

%       OutputPath: string that indicate where the data should be saved. If
%       not specified or left empty, the data will be saved in the folder
%       of the inputfile

%       NeuroBuffer: scalar indicating the duration of the time buffer in
%       ms that should be added before and after the onset and offset time
%       of vocalizations for extracting neural data. Default is the merged
%       threshold used in the extracted vocalization data

%       Flags: 2 scalar binary variable that indicate in the case of a CSC
%       input whether the Raw data and or the LFP should be extracted
%       Flags(1)= Raw data, Flags(2) = LFP in case a CSC

%       DenoiseT: binary variable, in case of a tetrode input file,
%       indicate whether the data should be extracted with various level of
%       denoising of the spikes (spikes automatically sorted by comparison
%       with Michael database of acceptable spikes)

%       Rthreshold: vector of values of correlation for the denoising of
%       spikes requested


% OUPUT
%       A file containing the neural data requested saved in OutputPath.
%       Spike arrival times are saved in ms in reference to the onset of
%       the vocalization bouts given by Voc_transc_time_refined.
if nargin<3
    NeuroBuffer=[];
end
if nargin<4 || isempty(Flags)
    Flags = [0 1];
end
if nargin<5
    DenoiseT = 0; % No sort of the tetrode spike from noise
end
if nargin<6
    Rthreshold = [0.92 0.94 0.96 0.98];
end
MaxEventDur = NaN; % Set to NaN: The neural data is extracted for the whole duration of each event
BaselineDur = 1000; % Duration of the baseline section in ms that is seeked at least BaselineDelay ms before sound onset
BaselineDelay = 5000;

%% Identify the neural data and run the corresponding extraction
[Path2Data, DataFile]=fileparts(InputDataFile);
PathParts = regexp(Path2Data, '/', 'split');
Loggers_dir = ['/' fullfile(PathParts{1:end-2})];
if nargin<2 || isempty(OutputPath)
    OutputPath = Loggers_dir;
end



if contains(DataFile, 'Tetrode')
    % this is a tetrode file
    % Get the date of the recording
    Idx_ = strfind(DataFile, '_');
    Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));
    
    % Get the tetrode ID
    NeuralInputID = DataFile(Idx_(end)+2);
    
    % Get the subject ID
    SubjectID = DataFile(1:5);
    
    % Loop through audio data
    AudioDir = dir(fullfile(Loggers_dir, sprintf('%s*VocExtractData.mat', Date(3:end)))); % These are all the results of vocalization localization for both operant conditioning and free session
    Nsession = length(AudioDir);
    ExpStartTimes = cell(Nsession,1);
    for nn=1:Nsession
        Idx_2 = strfind(AudioDir(nn).name, '_');
        ExpStartTimes{nn} = AudioDir(nn).name((Idx_2(1)+1) : (Idx_2(2)-1));
        load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date(3:end), ExpStartTimes{nn})), 'Voc_transc_time_refined');
        AudioDir2 = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date(3:end), ExpStartTimes{nn}))); % These are the results of vocalization identificaion and merging with Who calls for this session
        Idx_3 = strfind(AudioDir2(nn).name, '_');
        Idxmat = strfind(AudioDir2(nn).name, '.mat');
        if isempty(NeuroBuffer)
            NeuroBuffer = str2double(AudioDir2.name((Idx_3+1):(Idxmat-1))); % This is the merged threshold used in the extracted vocalization data in ms we want to use the same value as the neural buffer
        end
        
        % Make sure to re-order values in chronological order
        [~,SortInd] = sort(Voc_transc_time_refined(:,1));
        
        % Find the boundaries for obtaining silence sections of 1 second before 1s of each vocalization event
        BSL_transc_time_refined = find_dead_time(Voc_transc_time_refined(SortInd),BaselineDelay,BaselineDur);
        % extract Tetrode data
        [Voc_NeuroT] = extract_timeslot_Tetrode(InputDataFile, Voc_transc_time_refined(SortInd), BSL_transc_time_refined, NeuroBuffer,MaxEventDur, DenoiseT, Rthreshold, SubjectID, Date, NeuralInputID);
        OutputFile = fullfile(OutputPath, sprintf('%s_%s_%s_Tet%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID));
        if exist(OutputFile, 'file')
            save(OutputFile, 'Voc_NeuroT','-append');
        else
            save(OutputFile, 'Voc_NeuroT');
        end
    end
    
elseif contains(DataFile, 'CSC')
    % this is a Raw file
    % Get the date of the recording
    Idx_ = strfind(DataFile, '_');
    Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));
    
    % Get the channel ID
    NeuralInputID = DataFile((strfind(DataFile, 'CSC')+3):(strfind(DataFile, '.mat')-1));
    
    % Get the subject ID
    SubjectID = DataFile(1:5);
    
    % Loop through audio data
    AudioDir = dir(fullfile(Loggers_dir, sprintf('%s*VocExtractData.mat', Date(3:end)))); % These are all the results of vocalization localization for both operant conditioning and free session
    Nsession = length(AudioDir);
    ExpStartTimes = cell(Nsession,1);
    for nn=1:Nsession
        Idx_2 = strfind(AudioDir(nn).name, '_');
        ExpStartTimes{nn} = AudioDir(nn).name((Idx_2(1)+1) : (Idx_2(2)-1));
        load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date(3:end), ExpStartTimes{nn})), 'Voc_transc_time_refined');
        AudioDir2 = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date(3:end), ExpStartTimes{nn}))); % These are the results of vocalization identificaion and merging with Who calls for this session
        Idx_3 = strfind(AudioDir2(nn).name, '_');
        Idxmat = strfind(AudioDir2(nn).name, '.mat');
        if isempty(NeuroBuffer)
            NeuroBuffer = str2double(AudioDir2.name((Idx_3+1):(Idxmat-1))); % This is the merged threshold used in the extracted vocalization data in ms we want to use the same value as the neural buffer
        end
        
        % Make sure to re-order values in chronological order
        [~,SortInd] = sort(Voc_transc_time_refined(:,1));
        
        % Find the boundaries for obtaining silence sections of 1 second before 1s of each vocalization event
        BSL_transc_time_refined = find_dead_time(Voc_transc_time_refined(SortInd),BaselineDelay,BaselineDur);
        
        [Voc_NeuroRaw, Voc_NeuroLFP] = extract_timeslot_RawLFP(InputDataFile, Voc_transc_time_refined(SortInd), BSL_transc_time_refined, NeuroBuffer,MaxEventDur,NeuralInputID, Flags);
        if Flags(1)
            OutputFile = fullfile(OutputPath, sprintf('%s_%s_%s_Raw%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID));
            if exist(OutputFile, 'file')
                save(OutputFile, 'Voc_NeuroRaw','-append');
            else
                save(OutputFile, 'Voc_NeuroRaw');
            end
        end
        if Flags(2) % extract LFP data
            OutputFile = fullfile(OutputPath, sprintf('%s_%s_%s_LFP%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID));
            if exist(OutputFile, 'File')
                save(OutputFile, 'Voc_NeuroLFP', '-append');
            else
                save(OutputFile, 'Voc_NeuroLFP');
            end
        end
    end
elseif contains(DataFile,'SS')
    % this is a spike sorted unit
    % Get the date of the recording
    Idx_ = strfind(DataFile, '_');
    Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));
    AudioDataPath = ['/' fullfile(PathParts{1:end-4}, 'audio',Date)];
    
    % Get the tetrode ID
    NeuralInputID{1} = DataFile(strfind(DataFile, 'TT')+2);
    % Get the SS ID
    NeuralInputID{2} = DataFile((Idx_(end)+1):end);
    % Get the SS Q
    NeuralInputID{3} = DataFile(strfind(DataFile, '_SS')+3);
    
    % Get the subject ID
    SubjectID = DataFile(1:5);
    
    % Find if there is any period of unstability for the neural
    % activity %% NOT IMPLEMENTED AS OF NOW
    % load(fullfile(OutputPath, sprintf('%s_%s_SS%s_%s-%s.mat', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2}),'QualitySSU');
    
    % Loop through audio data
    AudioDir_all = dir(fullfile(Loggers_dir, sprintf('%s*VocExtractData*.mat', Date(3:end)))); % These are all the results of vocalization localization for both operant conditioning and free session
    TrueFiles = nan(length(AudioDir_all),1);
    FileOrder = nan(length(AudioDir_all),1);
    SessionTime_all = cell(length(AudioDir_all),1);
    for ff=1:length(AudioDir_all)
        TrueFiles(ff)=(length(strfind(AudioDir_all(ff).name, '_'))==2);
        SessionTime_all{ff} = AudioDir_all(ff).name(8:11);
        if TrueFiles(ff)
            Ind_ = strfind(AudioDir_all(ff).name, '.');
        else
            Ind_ = strfind(AudioDir_all(ff).name, '_');
            Ind_ = Ind_(end);
        end
        IndData = strfind(AudioDir_all(ff).name, 'Data');
        FileOrder(ff) = str2double(AudioDir_all(ff).name((IndData+4):(Ind_-1)));
    end
    AudioDir = AudioDir_all(logical(TrueFiles));
    FileOrder = FileOrder(logical(TrueFiles));
    SessionTime = SessionTime_all(logical(TrueFiles));
    ExpStartTimes = unique(SessionTime);
    Nsession = length(ExpStartTimes);
    for nn=1:Nsession
        % Retrieve the transctime of vocalizations in all files
        % corresponding to the session
        FileInd = find(contains(SessionTime, ExpStartTimes{nn}));
        [~,OrdInd] = sort(FileOrder(FileInd));
        AudioDir_local = AudioDir(FileInd);
        AudioDir_local = AudioDir_local(OrdInd);
        Voc_transc_time_refined_All = cell(1,length(FileInd));
        for ff=1:length(FileInd)
            load(fullfile(AudioDir_local(ff).folder,AudioDir_local(ff).name), 'Voc_transc_time_refined');
            Voc_transc_time_refined_All{ff} = Voc_transc_time_refined';
        end
        Voc_transc_time_refined = [Voc_transc_time_refined_All{:}]';
        
        load(fullfile(AudioDataPath, sprintf('%s_%s_VocExtractTimes.mat', Date(3:end), ExpStartTimes{nn})),'Re_transc_time');
        if ~exist('Re_transc_time','var')
            Re_transc_time = nan(size(Voc_transc_time_refined,1),1);
        end
        
        AudioDir2 = AudioDir_all(intersect(find(contains(SessionTime_all, ExpStartTimes{nn})), find(~TrueFiles))); % These are the results of vocalization identificaion and merging with Who calls for this session
        Idx_3 = strfind(AudioDir2(1).name, '_');
        Idxmat = strfind(AudioDir2(1).name, '.mat');
        if isempty(NeuroBuffer)
            NeuroBuffer = str2double(AudioDir2(1).name((Idx_3(end)+1):(Idxmat-1))); % This is the merged threshold used in the extracted vocalization data in ms we want to use the same value as the neural buffer
        end
        % Make sure to re-order values in chronological order
        [~,SortInd] = sort(Voc_transc_time_refined(:,1));
        % Find the boundaries for obtaining silence sections of BaselineDur second before BaselineDelay s of each vocalization event
        BSL_transc_time_refined = find_dead_time(Voc_transc_time_refined(SortInd,:),BaselineDelay,BaselineDur);
        % Extract Spike data
        [Voc_NeuroSSU] = extract_timeslot_SSU(InputDataFile, Voc_transc_time_refined(SortInd,:), BSL_transc_time_refined,Re_transc_time(SortInd), NeuroBuffer,MaxEventDur,Inf);
        Voc_NeuroSSU.BSL_transc_time_refined = BSL_transc_time_refined;
        Voc_NeuroSSU.SortInd=SortInd;
        
        OutputFile = fullfile(OutputPath, sprintf('%s_%s_%s_SS%s_%s-%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID{3},NeuralInputID{1},NeuralInputID{2}));
        if exist(OutputFile, 'file')
            save(OutputFile, 'Voc_NeuroSSU','-append');
        else
            save(OutputFile, 'Voc_NeuroSSU');
        end
        clear Re_transc_time Voc_transc_time_refined
    end
end


%% INTERNAL FUNCTIONS
% Calculate sections of times between behavioral events
% where we can extract baseline activity
    function [DeadTimes] = find_dead_time(InputTime, Delay, Duration)
        % InputTime is a 2 column vector with # rows = # behavioral events.
        % Delay: how long before the onset of each event in ms
        % Duration: how long should be the baseline extract in ms
        NEvent = size(InputTime,1);
        DeadTimes = nan(size(InputTime));
        [~, IndSort] = sort(InputTime(:,1));
        SortInputTime = InputTime(IndSort,:);
        for ee=1:NEvent
            PotentialOnset = SortInputTime(ee,1)-Duration-Delay;
            if ee==1 % Nothing to fear before
                DeadTimes(IndSort(ee),:) = PotentialOnset + [0 Duration];
            else
                if PotentialOnset>SortInputTime(ee-1,2) % all good to go
                    DeadTimes(IndSort(ee),:) = PotentialOnset + [0 Duration];
                elseif ((SortInputTime(ee,1)>SortInputTime(ee-1,1) && SortInputTime(ee,1)<SortInputTime(ee-1,2)) || (SortInputTime(ee,2)>SortInputTime(ee-1,1) && SortInputTime(ee,2)<SortInputTime(ee-1,2))) && ((PotentialOnset + Duration)<(SortInputTime(ee-1,1)+Delay)) % These extracts largely overlap, and the background is taken before the onset so good!  
                    DeadTimes(IndSort(ee),:) = PotentialOnset + [0 Duration];
                else%take at least Delay second after previous event offset and up to Delay second before current event onset
                    DeadTimes(IndSort(ee),1) = SortInputTime(ee-1,2) + Delay;
                    DeadTimes(IndSort(ee),2) = SortInputTime(ee,1) - Delay;
                    if diff(DeadTimes(IndSort(ee),:))<0 % No baseline available here!
                        DeadTimes(IndSort(ee),:) = nan(1,2);
                    end
                    
                end
            end
        end
    end



% Extracting tetrode data
    function [OutData] = extract_timeslot_Tetrode(InputFile, Voc_transc_time, BSL_transc_time, Buffer,MaxEventDur, DenoiseT, Rthreshold, SubjectID, Date, TID)
        Nevent = size(Voc_transc_time,1);
        OutData.SpikeTVoc = cell(Nevent,1);
        OutData.SpikeTVocDeNoiseInd = cell(Nevent,length(Rthreshold));
        OutData.SpikeTBSL = cell(Nevent,1);
        OutData.SpikeTBSLDeNoiseInd = cell(Nevent,length(Rthreshold));
        % reframe the extraction windows of each vocalization
        [EventOnset_time ,EventOffset_time] = reframe(Voc_transc_time, Buffer, MaxEventDur);
        % Load the spike arrival times for that tetrode
        [Path,~] = fileparts(InputFile);
        Spikes = load(fullfile(Path, sprintf('%s_%s_Tetrode_spikes_time_T%s.mat',SubjectID, Date, TID)), 'Spike_arrival_times');
        Snip = load(fullfile(Path, sprintf('%s_%s_Tetrode_spikes_snippets_T%s.mat',SubjectID, Date, TID)), 'Snippets');
        % loop through vocalizations and extract spike arrival times
        for vv=1:Nevent
            % Find the spike arrival times that are between the
            % requested times and center them to the onset of the
            % behavioral event
            SpikeT_local = logical((Spikes.Spike_arrival_times>(EventOnset_time(vv)*10^3)) .* (Spikes.Spike_arrival_times<(EventOffset_time(vv)*10^3)));
            OutData.SpikeTVoc{vv} = Spikes.Spike_arrival_times(SpikeT_local)/10^3 - Voc_transc_time(vv,1);
            if DenoiseT % get the indices of denoised spikes according to threshold(s)
                for rr=1:length(Rthreshold)
                    OutData.SpikeTVocDeNoiseInd{vv,rr} = sort_spike_from_noise(Snip.Snippets(:,:,SpikeT_local), Rthreshold(rr));
                end
            end
            SpikeT_local = logical((Spikes.Spike_arrival_times>(BSL_transc_time(vv,1)*10^3)) .* (Spikes.Spike_arrival_times<(BSL_transc_time(vv,2)*10^3)));
            OutData.SpikeTBSL{vv} = Spikes.Spike_arrival_times(SpikeT_local)/10^3 - Voc_transc_time(vv,1);
            if DenoiseT % get the indices of denoised spikes according to threshold(s)
                for rr=1:length(Rthreshold)
                    OutData.SpikeTBSLDeNoiseInd{vv,rr} = sort_spike_from_noise(Snip.Snippets(:,:,SpikeT_local), Rthreshold(rr));
                end
            end
        end
    end





% Extracting Raw and/or LFP data
    function [OutDataRaw, OutDataLFP] = extract_timeslot_RawLFP(InputFile, Voc_transc_time, BSL_transc_time, Buffer,MaxEventDur,ChannelID, Flag)
        Nevent = size(Voc_transc_time,1);
        OutDataRaw = struct();
        OutDataLFP = struct();
        if Flag(1)
            OutDataRaw.RawVoc = cell(Nevent,1);
            OutDataRaw.RawBSL = cell(Nevent,1);
        end
        if Flag(2)
            OutDataLFP.LFPVoc = cell(Nevent,1);
            OutDataLFP.LFPBSL = cell(Nevent,1);
        end
        % reframe the extraction windows of each vocalization
        [EventOnset_time ,EventOffset_time] = reframe(Voc_transc_time, Buffer, MaxEventDur);
        % load the data
        LData = load(InputFile, 'Timestamps_of_first_samples_usec', 'Estimated_channelFS_Transceiver','AD_count_int16','Indices_of_first_and_last_samples','Sampling_period_usec_Logger');
        % loop through events and extract the snippet of Raw data
        NanFSInd = find(isnan(LData.Estimated_channelFS_Transceiver));
        GoodFSInd = find(~isnan(LData.Estimated_channelFS_Transceiver));
        for vv=1:Nevent
            fprintf('Raw data Channel %d Event %d/%d\n', ChannelID, vv, Nevent);
            if prod(~isnan(Voc_transc_time(vv,:)))
                % find the time stamp on the logger that is closest to before
                % the snippet of event onset
                [IndSampOn, IndSampOff, FS] = time2indices_logger(EventOnset_time(vv), EventOffset_time(vv), length(LData.AD_count_int16), LData.Timestamps_of_first_samples_usec, LData.Indices_of_first_and_last_samples, LData.Estimated_channelFS_Transceiver, LData.Sampling_period_usec_Logger, NanFSInd, GoodFSInd);
                % extract the voltage snippet and bandpass the raw signal
                % if extracting LFP
                if ~isempty(IndSampOn) && Flag(1) && ~Flag(2)
                    OutDataRaw.RawVoc{vv} = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                elseif ~isempty(IndSampOn) && ~Flag(1) && Flag(2)
                    [z,p,k] = butter(12,BandPassFilter(1)/(FS/2),'low');
                    sos = zp2sos(z,p,k);
                    RawVoc = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                    OutDataLFP.LFPVoc{vv} = (filtfilt(sos,1,RawVoc)); % % low-pass filter the voltage trace
                elseif ~isempty(IndSampOn) && Flag(1) && Flag(2)
                    OutDataRaw.RawVoc{vv} = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                    [z,p,k] = butter(12,BandPassFilter(1)/(FS/2),'low');
                    sos = zp2sos(z,p,k);
                    OutDataLFP.LFPVoc{vv} = (filtfilt(sos,1,OutDataRaw.RawVoc{vv})); % % low-pass filter the voltage trace
                end
                
            end
            
            if prod(~isnan(BSL_transc_time(vv,:)))
                % find the time stamp on the logger that is closest to before
                % the snippet of event onset
                [IndSampOn, IndSampOff, FS] = time2indices_logger(BSL_transc_time(vv,1), BSL_transc_time(vv,2), length(LData.AD_count_int16), LData.Timestamps_of_first_samples_usec, LData.Indices_of_first_and_last_samples, LData.Estimated_channelFS_Transceiver, LData.Sampling_period_usec_Logger, NanFSInd, GoodFSInd);
                % extract the voltage snippet and bandpass the raw signal
                % if extracting LFP
                if ~isempty(IndSampOn) && Flag(1) && ~Flag(2)
                    OutDataRaw.RawBSL{vv} = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                elseif ~isempty(IndSampOn) && ~Flag(1) && Flag(2)
                    [z,p,k] = butter(12,BandPassFilter(1)/(FS/2),'low');
                    sos = zp2sos(z,p,k);
                    RawBSL = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                    OutDataLFP.LFPBSL{vv} = (filtfilt(sos,1,RawBSL)); % % low-pass filter the voltage trace
                elseif ~isempty(IndSampOn) && Flag(1) && Flag(2)
                    OutDataRaw.RawBSL{vv} = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                    [z,p,k] = butter(12,BandPassFilter(1)/(FS/2),'low');
                    sos = zp2sos(z,p,k);
                    OutDataLFP.LFPBSL{vv} = (filtfilt(sos,1,OutDataRaw.RawBSL{vv})); % % low-pass filter the voltage trace
                end
                
            end
        end
    end



% Extracting spike sorted unit data
    function [OutData] = extract_timeslot_SSU(InputFile, Voc_transc_time, BSL_transc_time,Re_transc_time, Buffer,MaxEventDur,DeadTime_usec)
        % Don't extract data for events that happened after the DeadTime of
        % the unit
        VocIdxSSUAlive = find(sum(Voc_transc_time<(DeadTime_usec.*10^-3),2)==2);
        Nevent = length(VocIdxSSUAlive);
        OutData.VocIdxSSUAlive = VocIdxSSUAlive;
        OutData.SpikeSUVoc = cell(Nevent,1);
        OutData.SpikeSUBSL = cell(Nevent,1);
        OutData.ReTime = nan(Nevent,1);
        % reframe the extraction windows of each vocalization
        [EventOnset_time ,EventOffset_time] = reframe(Voc_transc_time, Buffer, MaxEventDur);
        % loading the single unit spike arrival times
        Spikes = load(InputFile, 'Spike_arrival_times');
        % loop through vocalizations and extract spike arrival times
        for vv=1:Nevent
            ii = VocIdxSSUAlive(vv);
            % Find the spike arrival times that are between the
            % requested times and center them to the onset of the
            % behavioral event, save in ms
            OutData.SpikeSUVoc{vv} = Spikes.Spike_arrival_times(logical((Spikes.Spike_arrival_times>(EventOnset_time(ii)*10^3)) .* (Spikes.Spike_arrival_times<(EventOffset_time(ii)*10^3))))/10^3 - Voc_transc_time(ii,1);
            OutData.SpikeSUBSL{vv} = Spikes.Spike_arrival_times(logical((Spikes.Spike_arrival_times>(BSL_transc_time(ii,1)*10^3)) .* (Spikes.Spike_arrival_times<(BSL_transc_time(ii,2)*10^3))))/10^3 - Voc_transc_time(ii,1);
            if sum(isnan(Re_transc_time)) == length(Re_transc_time) % there is no reward time
                OutData.ReTime(vv) = NaN;
            else
                if vv<Nevent
                    vv_timecheck = find(logical((Re_transc_time>Voc_transc_time(ii,1)) .* (Re_transc_time<Voc_transc_time(ii+1,1))));
                else
                    vv_timecheck = find(Re_transc_time>(Voc_transc_time(ii,1)));
                end
                if any(vv == vv_timecheck)
                    OutData.ReTime(vv) = Re_transc_time(vv) - Voc_transc_time(ii,1);
                elseif isinf(Re_transc_time(vv))
                    OutData.ReTime(vv) = Re_transc_time(vv) - Voc_transc_time(ii,1);
                else
                    warning('Discrepancy between the guessed reward time given the vocalization onset times and the order in which they are in Re_transc_time\n')
                    keyboard
                    GuessVV = Re_transc_time(vv) - Voc_transc_time(ii,1);
                    GuessVVTC = Re_transc_time(vv_timecheck) - Voc_transc_time(ii,1);
                    if isnan(GuessVV) && isempty(GuessVVTC) || (GuessVV)<0 && isempty(GuessVVTC)
                        OutData.ReTime(vv) = NaN; % No reward for that vocalization
                    elseif  isempty(GuessVVTC)
                        OutData.ReTime(vv) = NaN; % No reward for that vocalization
                    elseif GuessVV<0 && GuessVVTC>0
                        OutData.ReTime(vv) = GuessVVTC;
                    elseif GuessVV>0 && GuessVVTC<0
                        OutData.ReTime(vv) = GuessVV;
                    else
                        keyboard
                    end
                end
                
            end
        end
    end



% Enlarge the window of extraction
    function [EventOnset_time, EventOffset_time] = reframe(OnsetOffset_time, Buffer, MaxEventDur)
        % reframe in time adding Buffer ms before and after the
        % OnsetOffset_time given in ms to the extent of creating an event
        % of maximum duration MaxEventDur. No maximum duration implemented
        % if MaxEventDur is Nan.
        Nevent = size(OnsetOffset_time,1);
        EventOnset_time = nan(Nevent,1);
        EventOffset_time = nan(Nevent, 1);
        for vv=1:Nevent
            if prod(~isnan(OnsetOffset_time(vv,:)))
                EventOnset_time(vv) = OnsetOffset_time(vv,1) - Buffer;
                if isnan(MaxEventDur)
                    EventOffset_time(vv) = OnsetOffset_time(vv,2) + Buffer;
                elseif diff(OnsetOffset_time(vv,:))>(MaxEventDur + Buffer)
                    EventOffset_time(vv) = EventOnset_time(vv) + MaxEventDur + 2*Buffer;
                elseif diff(OnsetOffset_time(vv,:))<MaxEventDur
                    EventOffset_time(vv) = OnsetOffset_time(vv,2) + Buffer;
                elseif diff(OnsetOffset_time(vv,:))>MaxEventDur
                    EventOffset_time(vv) = EventOnset_time(vv) + MaxEventDur + 2*Buffer;
                end
            end
        end
    end



% Find the indices onset/offset in the raw logger file corresponding to requested transceiver times 
    function [IndSampOn, IndSampOff, FS] = time2indices_logger(TimeOnset, TimeOffset, RawDataLength, Timestamps_of_first_samples_usec, Indices_of_first_and_last_samples,Estimated_channelFS_Transceiver,Sampling_period_usec_Logger,NanFSInd, GoodFSInd)
        IndTSOn = find(Timestamps_of_first_samples_usec<(TimeOnset*10^3), 1, 'Last');
        if ~isempty(IndTSOn) %This event did happen after recording onset
            % deduct the corresponding onset  sample
            if ~isempty(intersect(NanFSInd,IndTSOn)) || IndTSOn>length(Estimated_channelFS_Transceiver) % there is no sample frequency estimate for that recorded file, take the previous or following estimate
                OKFs = [find(GoodFSInd<IndTSOn,1, 'Last') find(GoodFSInd>IndTSOn,1, 'First')];
                if ~isempty(OKFs)
                    Local_FS = Estimated_channelFS_Transceiver(GoodFSInd(OKFs(1)));
                else
                    Local_FS = 1/Sampling_period_usec_Logger*10^6;
                end
            else
                Local_FS =Estimated_channelFS_Transceiver(IndTSOn);
            end
            IndSampOn = round(Indices_of_first_and_last_samples(IndTSOn,1) + Local_FS*(10^-6)*(TimeOnset*10^3 - Timestamps_of_first_samples_usec(IndTSOn)));
            
            % find the time stamp on the logger that is closest to after
            % the snippet of event offset
            if (TimeOffset*10^3)<Timestamps_of_first_samples_usec(end) % The end of the event is before the onset of the last recorded file
                IndTSOff = find(Timestamps_of_first_samples_usec>(TimeOffset*10^3), 1, 'First');
                % deduct the corresponding onset offset samples
                if ~isempty(intersect(NanFSInd,IndTSOff)) || IndTSOff>length(Estimated_channelFS_Transceiver) % there is no sample frequency estimate for that recorded file, take the previous or following estimate
                    OKFs = [find(GoodFSInd<IndTSOff,1, 'Last') find(GoodFSInd>IndTSOff,1, 'First')];
                    if ~isempty(OKFs)
                        Local_FS = Estimated_channelFS_Transceiver(GoodFSInd(OKFs(1)));
                        FS = nanmean(Estimated_channelFS_Transceiver(IndTSOn:GoodFSInd(OKFs(1))));
                    else
                        Local_FS = 1/Sampling_period_usec_Logger*10^6;
                    end
                    
                else
                    Local_FS = Estimated_channelFS_Transceiver(IndTSOff);
                    FS = nanmean(Estimated_channelFS_Transceiver(IndTSOn:IndTSOff));
                end
                IndSampOff = round(Indices_of_first_and_last_samples(IndTSOff,1) - Local_FS*(10^-6)*(Timestamps_of_first_samples_usec(IndTSOff) - EventOffset_time(vv)*10^3));
                
            else % The end of the event is after the onset of the last recorded file
                IndTSOff = length(Timestamps_of_first_samples_usec);
                % deduct the corresponding offset samples
                if ~isempty(intersect(NanFSInd,IndTSOff)) || IndTSOff>length(Estimated_channelFS_Transceiver) % there is no sample frequency estimate for that recorded file, take the previous or following estimate
                    OKFs = [find(GoodFSInd<IndTSOff,1, 'Last') find(GoodFSInd>IndTSOff,1, 'First')];
                    if ~isempty(OKFs)
                        Local_FS = Estimated_channelFS_Transceiver(GoodFSInd(OKFs(1)));
                    else
                        Local_FS = 1/Sampling_period_usec_Logger*10^6;
                    end
                else
                    Local_FS =Estimated_channelFS_Transceiver(IndTSOff);
                end
                IndSampOff = round(Indices_of_first_and_last_samples(IndTSOff,1) + Local_FS*(10^-6)*(TimeOffset*10^3 - Timestamps_of_first_samples_usec(IndTSOff)));
                FS = nanmean(Estimated_channelFS_Transceiver(IndTSOn:end));
            end
            IndSampOff = min(RawDataLength, IndSampOff);
        else
            IndSampOn = [];
            IndSampOff = [];
        end
    end

end
