function cut_neuralData_behav_perfile(InputDataFile, OutputPath, BehaviorType,Flags, DenoiseT, Rthreshold)
%% This function uses the calculation of behavior onset/offset in transceiver time (ms) calculated by get_logger_data_behav
% (AllActions_Time) And extract the corresponding neural data in
% the input datafile.
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

%       BehaviorType: cell array of strings correcponding to the behaviors
%       that should be extracted ('licking' 'chewing' 'quiet' 'teeth-cleaning')

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
%       the behaviors given by AllActions_Time.

if nargin < 3
    BehaviorType = {'licking' 'chewing' 'quiet' 'teeth-cleaning'};
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
    
    % Loop through behavioral data files
    AudioDir = dir(fullfile(Loggers_dir, sprintf('%s*BehavExtractData.mat', Date(3:end)))); % These are all the results of behavior localization for free session
    Nsession = length(AudioDir);
    ExpStartTimes = cell(Nsession,1);
    for nn=1:Nsession
        Idx_2 = strfind(AudioDir(nn).name, '_');
        ExpStartTimes{nn} = AudioDir(nn).name((Idx_2(1)+1) : (Idx_2(2)-1));
        load(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData.mat', Date(3:end), ExpStartTimes{nn})), 'AllActions_Time', 'AllActions_ID','UActionText');
        % Select behavioral data
        [Behav_transc_time, Behav_What,Behav_Who] = sort_behavior(AllActions_Time, AllActions_ID, UActionText,BehaviorType);
        % extract Tetrode data
        [Behav_NeuroT] = extract_timeslot_Tetrode(InputDataFile, Behav_transc_time, DenoiseT, Rthreshold, SubjectID, Date, NeuralInputID);
        OutputFile = fullfile(OutputPath, sprintf('%s_%s_%s_Tet%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID));
        if exist(OutputFile, 'file')
            save(OutputFile, 'Behav_NeuroT','Behav_What','Behav_Who','-append');
        else
            save(OutputFile, 'Behav_NeuroT','Behav_What','Behav_Who');
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
    
    % Loop through behavioral data files
    AudioDir = dir(fullfile(Loggers_dir, sprintf('%s*BehavExtractData.mat', Date(3:end)))); % These are all the results of behavior localization for free session
    Nsession = length(AudioDir);
    ExpStartTimes = cell(Nsession,1);
    for nn=1:Nsession
        Idx_2 = strfind(AudioDir(nn).name, '_');
        ExpStartTimes{nn} = AudioDir(nn).name((Idx_2(1)+1) : (Idx_2(2)-1));
        load(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData.mat', Date(3:end), ExpStartTimes{nn})), 'AllActions_Time', 'AllActions_ID','UActionText');
        % Select behavioral data
        [Behav_transc_time, Behav_What,Behav_Who] = sort_behavior(AllActions_Time, AllActions_ID, UActionText,BehaviorType);
        % extract the channel data
        [Behav_NeuroRaw, Behav_NeuroLFP] = extract_timeslot_RawLFP(InputDataFile, Behav_transc_time, NeuralInputID, Flags);
        if Flags(1)
            OutputFile = fullfile(OutputPath, sprintf('%s_%s_%s_Raw%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID));
            if exist(OutputFile, 'file')
                save(OutputFile, 'Behav_NeuroRaw','Behav_What','Behav_Who','-append');
            else
                save(OutputFile, 'Behav_NeuroRaw', 'Behav_What','Behav_Who');
            end
        end
        if Flags(2) % extract LFP data
            OutputFile = fullfile(OutputPath, sprintf('%s_%s_%s_LFP%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID));
            if exist(OutputFile, 'File')
                save(OutputFile, 'Behav_NeuroLFP','Behav_What','Behav_Who', '-append');
            else
                save(OutputFile, 'Behav_NeuroLFP','Behav_What','Behav_Who');
            end
        end
    end
elseif contains(DataFile,'SS')
    % this is a spike sorted unit
    % Get the date of the recording
    Idx_ = strfind(DataFile, '_');
    Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));
    
    % Get the tetrode ID
    NeuralInputID{1} = DataFile(strfind(DataFile, 'TT')+2);
    % Get the SS ID
    NeuralInputID{2} = DataFile((Idx_(end)+1):end);
    
    % Get the subject ID
    SubjectID = DataFile(1:5);
    
    % Find if there is any period of unstability for the neural
    % activity %% THIS IS NOT IMPLEMENTED AS OF NOW, SHOULD BE AN INPUT OF
    % extract_timeslot_SSU
%     load(fullfile(OutputPath, sprintf('%s_%s_SSU%s-%s.mat', SubjectID, Date,NeuralInputID{1},NeuralInputID{2})),'QualitySSU');
    
    % Loop through behavioral data files
    AudioDir = dir(fullfile(Loggers_dir, sprintf('%s*BehavExtractData.mat', Date(3:end)))); % These are all the results of non-vocal Behavior localization for free session
    Nsession = length(AudioDir);
    ExpStartTimes = cell(Nsession,1);
    for nn=1:Nsession
        Idx_2 = strfind(AudioDir(nn).name, '_');
        ExpStartTimes{nn} = AudioDir(nn).name((Idx_2(1)+1) : (Idx_2(2)-1));
        load(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData.mat', Date(3:end), ExpStartTimes{nn})), 'AllActions_Time', 'AllActions_ID','UActionText');
        % Select behavioral data
        [Behav_transc_time, Behav_What,Behav_Who] = sort_behavior(AllActions_Time, AllActions_ID, UActionText,BehaviorType);
        % extract single unit spikes
        [Behav_NeuroSSU] = extract_timeslot_SSU(InputDataFile, Behav_transc_time,Inf);
        OutputFile = fullfile(OutputPath, sprintf('%s_%s_%s_SSU%s-%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID{1},NeuralInputID{2}));
        if exist(OutputFile, 'file')
            save(OutputFile, 'Behav_NeuroSSU','Behav_What','Behav_Who','-append');
        else
            save(OutputFile, 'Behav_NeuroSSU','Behav_What','Behav_Who');
        end
    end
end


%% INTERNAL FUNCTIONS

% Organize onset/offset time of behaviors as a 2 column matrix, indicate
% the type of behavior (Behav_what) and the identity of the bat performing the
% behavior
    function [Behav_transc_time, Behav_What,Behav_Who] = sort_behavior(AllActions_Time, AllActions_ID, UActionText,BehaviorType)
        % Get each behavior number of cuts to allocate space
        Num_Slots = nan(size(BehaviorType));
        IB_Ind = nan(size(BehaviorType));
        for bb=1:length(BehaviorType)
            IB_Ind_logical = contains(UActionText, BehaviorType{bb});
            if any(IB_Ind_logical)
                IB_Ind(bb) = find(IB_Ind_logical);
                Num_Slots(bb) = size(AllActions_Time{IB_Ind(bb)},1);
            end
        end
        % Initialize variables
        Behav_transc_time = nan(nansum(Num_Slots),2);
        Behav_What = cell(nansum(Num_Slots),1);
        Behav_Who =  nan(nansum(Num_Slots),1);
        BehavCount = 1;
        for bb=1:length(BehaviorType)
            if ~isnan(Num_Slots(bb))
                Behav_transc_time(BehavCount:(BehavCount+Num_Slots(bb)-1),:) = AllActions_Time{IB_Ind(bb)};
                Behav_What(BehavCount:(BehavCount+Num_Slots(bb)-1)) = BehaviorType(bb);
                Behav_Who(BehavCount:(BehavCount+Num_Slots(bb)-1),:) = AllActions_ID{IB_Ind(bb)};
                BehavCount = BehavCount + Num_Slots(bb);
            end
        end
    end



% Extracting tetrode data
    function [OutData] = extract_timeslot_Tetrode(InputFile, Voc_transc_time, DenoiseT, Rthreshold, SubjectID, Date, TID)
        Nevent = size(Voc_transc_time,1);
        OutData.SpikeTBehav = cell(Nevent,1);
        OutData.SpikeTBehavDeNoiseInd = cell(Nevent,length(Rthreshold));
        OutData.Duration = nan(Nevent,1);
        % Load the spike arrival times for that tetrode
        [Path,~] = fileparts(InputFile);
        Spikes = load(fullfile(Path, sprintf('%s_%s_Tetrode_spikes_time_T%s.mat',SubjectID, Date, TID)), 'Spike_arrival_times');
        Snip = load(fullfile(Path, sprintf('%s_%s_Tetrode_spikes_snippets_T%s.mat',SubjectID, Date, TID)), 'Snippets');
        % loop through vocalizations and extract spike arrival times
        for vv=1:Nevent
            % Find the spike arrival times that are between the
            % requested times and center them to the onset of the
            % behavioral event
            SpikeT_local = logical((Spikes.Spike_arrival_times>(Voc_transc_time(vv,1)*10^3)) .* (Spikes.Spike_arrival_times<(Voc_transc_time(vv,2)*10^3)));
            OutData.SpikeTBehav{vv} = Spikes.Spike_arrival_times(SpikeT_local)/10^3 - Voc_transc_time(vv,1);
            OutData.Duration(vv) = (Voc_transc_time(vv,2)*10^3)-(Voc_transc_time(vv,1)*10^3);
            if DenoiseT % get the indices of denoised spikes according to threshold(s)
                for rr=1:length(Rthreshold)
                    OutData.SpikeTBehavDeNoiseInd{vv,rr} = sort_spike_from_noise(Snip.Snippets(:,:,SpikeT_local), Rthreshold(rr));
                end
            end
        end
    end





% Extracting Raw and/or LFP data
    function [OutDataRaw, OutDataLFP] = extract_timeslot_RawLFP(InputFile, Voc_transc_time, ChannelID, Flag)
        Nevent = size(Voc_transc_time,1);
        OutDataRaw = struct();
        OutDataLFP = struct();
        if Flag(1)
            OutDataRaw.RawBehav = cell(Nevent,1);
        end
        if Flag(2)
            OutDataLFP.LFPBehav = cell(Nevent,1);
        end
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
                [IndSampOn, IndSampOff, FS] = time2indices_logger(Voc_transc_time(vv,1), Voc_transc_time(vv,2), length(LData.AD_count_int16), LData.Timestamps_of_first_samples_usec, LData.Indices_of_first_and_last_samples, LData.Estimated_channelFS_Transceiver, LData.Sampling_period_usec_Logger, NanFSInd, GoodFSInd);
                % extract the voltage snippet and bandpass the raw signal
                % if extracting LFP
                if ~isempty(IndSampOn) && Flag(1) && ~Flag(2)
                    OutDataRaw.RawBehav{vv} = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                elseif ~isempty(IndSampOn) && ~Flag(1) && Flag(2)
                    [z,p,k] = butter(12,BandPassFilter(1)/(FS/2),'low');
                    sos = zp2sos(z,p,k);
                    RawVoc = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                    OutDataLFP.LFPBehav{vv} = (filtfilt(sos,1,RawVoc)); % % low-pass filter the voltage trace
                elseif ~isempty(IndSampOn) && Flag(1) && Flag(2)
                    OutDataRaw.RawBehav{vv} = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                    [z,p,k] = butter(12,BandPassFilter(1)/(FS/2),'low');
                    sos = zp2sos(z,p,k);
                    OutDataLFP.LFPBehav{vv} = (filtfilt(sos,1,OutDataRaw.RawBehav{vv})); % % low-pass filter the voltage trace
                end
                
            end
        end
    end



% Extracting spike sorted unit data
    function [OutData] = extract_timeslot_SSU(InputFile, Behav_transc_time,DeadTime_usec)
        % Don't extract data for events that happened after the DeadTime of
        % the unit
        BehavIdxSSUAlive = find(sum(Behav_transc_time<(DeadTime_usec.*10^-3),2)==2);
        Nevent = length(BehavIdxSSUAlive);
        OutData.BehavIdxSSUAlive = BehavIdxSSUAlive;
        OutData.SpikeSUBehav = cell(Nevent,1);
        OutData.Duration = nan(Nevent,1);
        
        % loading the single unit spike arrival times
        Spikes = load(InputFile, 'Spike_arrival_times');
        % loop through vocalizations and extract spike arrival times
        for vv=1:Nevent
            ii = BehavIdxSSUAlive(vv);
            % Find the spike arrival times that are between the
            % requested times and center them to the onset of the
            % behavioral event, save in ms
            OutData.SpikeSUBehav{vv} = Spikes.Spike_arrival_times(logical((Spikes.Spike_arrival_times>(Behav_transc_time(ii,1)*10^3)) .* (Spikes.Spike_arrival_times<(Behav_transc_time(ii,2)*10^3))))/10^3 - Behav_transc_time(ii,1);
            OutData.Duration(vv) = Behav_transc_time(ii,2)*10^3 - Behav_transc_time(ii,1)*10^3; % Duration value in ms
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
