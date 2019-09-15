function cut_neuralData_voc_perfile(InputDataFile,  ExpStartTime, Flags, DenoiseT, Rthreshold)
%% This function uses the better estimation of vocalization onset/offset in transceiver time (ms) calculated by get_logger_data_voc
% (Voc_transc_time_refined) And extract the corresponding neural data in
% the inout datafile as long as neural data of baseline activity in the
% seconds preceding the vocalization

% INPUT:
%       InputDataFile: path to a tetrode mat file, a single unit mat file or the
%       raw data (*CSC*.mat file)

%       Flags: 2 scalar binary variable that indicate in the case of a CSC
%       input whether the Raw data and or the LFP should be extracted
%       Flags(1)= Raw data, Flags(2) = LFP in case a CSC

%       DenoiseT: binary variable, in case of a tetrode input file,
%       indicate whether the data should be extracted with various level of
%       denoising of the spikes (spikes automatically sorted by comparison
%       with Michael database of acceptable spikes)

%       Rthreshold: vector of values of correlation for the denoising of
%       spikes requested


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
BaselineDur = 1000; % Duration of the baseline section in ms that is seeked at least BaselineDelay ms before sound onset
BaselineDelay = 1000;

%% Identify the neural data and run the corresponding extraction
[Path2Data, DataFile]=fileparts(InputDataFile);
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
    PathParts = regexp(Path2Data, '/', 'split');
    Loggers_dir = fullfile(PathParts(1:end-2));
    AudioDir = dir(fullfile(Loggers_dir, sprintf('%s*VocExtractData.mat', Date))); % These are all the results of vocalization localization for both operant conditioning and free session
    Nsession = length(AudioDir);
    ExpStartTimes = cell(Nsession,1);
    for nn=1:Nsession
        Idx_2 = strfind(AudioDir(nn).name, '_');
        ExpStartTimes{nn} = AudioDir(nn).name((Idx_2(1)+1) : (Idx_2(2)-1));
        load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTimes{nn})), 'Voc_transc_time_refined');
        AudioDir2 = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date, ExpStartTimes{nn}))); % These are the results of vocalization identificaion and merging with Who calls for this session
        Idx_3 = strfind(AudioDir2(nn).name, '_');
        Idxmat = strfind(AudioDir2(nn).name, '.mat');
        NeuroBuffer = str2double(AudioDir2.name((Idx_3+1):(Idxmat-1))); % This is the merged threshold used in the extracted vocalization data in ms we want to use the same value as the neural buffer
        % Find the boundaries for obtaining silence sections of 1 second before 1s of each vocalization event
        BSL_transc_time_refined = find_dead_time(Voc_transc_time_refined,BaselineDelay,BaselineDur);
        % extract Tetrode data
        [Voc_NeuroT] = extract_timeslot_Tetrode(InputDataFile, Voc_transc_time_refined, BSL_transc_time_refined, NeuroBuffer,MaxEventDur, DenoiseT, Rthreshold, SubjectID, Date, NeuralInputID);
        save(fullfile(Loggers_dir, sprintf('%s_%s_%s_Tet%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID)), Voc_NeuroT);
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
    PathParts = regexp(Path2Data, '/', 'split');
    Loggers_dir = fullfile(PathParts(1:end-2));
    AudioDir = dir(fullfile(Loggers_dir, sprintf('%s*VocExtractData.mat', Date))); % These are all the results of vocalization localization for both operant conditioning and free session
    Nsession = length(AudioDir);
    ExpStartTimes = cell(Nsession,1);
    for nn=1:Nsession
        Idx_2 = strfind(AudioDir(nn).name, '_');
        ExpStartTimes{nn} = AudioDir(nn).name((Idx_2(1)+1) : (Idx_2(2)-1));
        load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTimes{nn})), 'Voc_transc_time_refined');
        AudioDir2 = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date, ExpStartTimes{nn}))); % These are the results of vocalization identificaion and merging with Who calls for this session
        Idx_3 = strfind(AudioDir2(nn).name, '_');
        Idxmat = strfind(AudioDir2(nn).name, '.mat');
        NeuroBuffer = str2double(AudioDir2.name((Idx_3+1):(Idxmat-1))); % This is the merged threshold used in the extracted vocalization data in ms we want to use the same value as the neural buffer
        % Find the boundaries for obtaining silence sections of 1 second before 1s of each vocalization event
        BSL_transc_time_refined = find_dead_time(Voc_transc_time_refined,BaselineDelay,BaselineDur);
        
        if Flags(1) % extract Raw data
            [Voc_NeuroRaw] = extract_timeslot_Raw(InputDataFile, Voc_transc_time_refined, BSL_transc_time_refined, NeuroBuffer,MaxEventDur);
            save(fullfile(Loggers_dir, sprintf('%s_%s_%s_Raw%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID)), Voc_NeuroRaw);
        end
        if Flags(2) % extract LFP data
            [Voc_NeuroLFP] = extract_timeslot_LFP(InputDataFile, Voc_transc_time_refined, BSL_transc_time_refined, NeuroBuffer,MaxEventDur);
            save(fullfile(Loggers_dir, sprintf('%s_%s_%s_LFP%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID)), Voc_NeuroLFP);
        end
    end
elseif contains('SS')
    % this is a spike sorted unit
    % Get the date of the recording
    Idx_ = strfind(DataFile, '_');
    Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));
    
    % Get the tetrode ID
    NeuralInputID{1} = DataFile(strfind(DataFile, 'TT')+3);
    % Get the SS ID
    NeuralInputID{2} = DataFile((strfind(DataFile, 'SS')+2):(strfind(DataFile, '.mat')-1));
    
    % Get the subject ID
    SubjectID = DataFile(1:5);
    
    % Loop through audio data
    PathParts = regexp(Path2Data, '/', 'split');
    Loggers_dir = fullfile(PathParts(1:end-2));
    AudioDir = dir(fullfile(Loggers_dir, sprintf('%s*VocExtractData.mat', Date))); % These are all the results of vocalization localization for both operant conditioning and free session
    Nsession = length(AudioDir);
    ExpStartTimes = cell(Nsession,1);
    for nn=1:Nsession
        Idx_2 = strfind(AudioDir(nn).name, '_');
        ExpStartTimes{nn} = AudioDir(nn).name((Idx_2(1)+1) : (Idx_2(2)-1));
        load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTimes{nn})), 'Voc_transc_time_refined');
        AudioDir2 = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date, ExpStartTimes{nn}))); % These are the results of vocalization identificaion and merging with Who calls for this session
        Idx_3 = strfind(AudioDir2(nn).name, '_');
        Idxmat = strfind(AudioDir2(nn).name, '.mat');
        NeuroBuffer = str2double(AudioDir2.name((Idx_3+1):(Idxmat-1))); % This is the merged threshold used in the extracted vocalization data in ms we want to use the same value as the neural buffer
        % Find the boundaries for obtaining silence sections of 1 second before 1s of each vocalization event
        BSL_transc_time_refined = find_dead_time(Voc_transc_time_refined,BaselineDelay,BaselineDur);
        [Voc_NeuroSSU] = extract_timeslot_SSU(InputDataFile, Voc_transc_time_refined, BSL_transc_time_refined, NeuroBuffer,MaxEventDur);
        save(fullfile(Loggers_dir, sprintf('%s_%s_%s_SSU%s-%s.mat', SubjectID, Date, ExpStartTimes{nn},NeuralInputID{1},NeuralInputID{2})), Voc_NeuroSSU);
    end
end


%% INTERNAL FUNCTIONS
% Internal function to calculate sections of times where we can extract baseline activity
    function [DeadTimes] = find_dead_time(InputTime, Delay, Duration)
        NEvent = size(InputTime,1);
        DeadTimes = nan(size(InputTime));
        for ee=1:NEvent
            PotentialOnset = InputTime(ee,1)-Duration-Delay;
            if ee==1 % Nothing to fear before
                DeadTimes(ee,:) = PotentialOnset + [0 Duration];
            else
                if PotentialOnset>InputTime(ee-1,2) % all good to go
                    DeadTimes(ee,:) = PotentialOnset + [0 Duration];
                else %take at least one second after previous event offset and up to one second before current event onset
                    DeadTimes(ee,1) = InputTime(ee-1,2) + Delay;
                    DeadTimes(ee,2) = InputTime(ee,1) - Delay;
                    if diff(DeadTimes(ee,:))<0 % No baseline available here!
                        DeadTimes(ee,:) = nan(1,2);
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





% Extracting Raw data
    function [OutData] = extract_timeslot_Raw(InputFile, Voc_transc_time, BSL_transc_time, Buffer,MaxEventDur)
        
    end



% Extracting LFP data
    function [OutData] = extract_timeslot_LFP(InputFile, Voc_transc_time, BSL_transc_time, Buffer,MaxEventDur)
        
    end



% Extracting spike sorted unit data
    function [OutData] = extract_timeslot_SSU(InputFile, Voc_transc_time, BSL_transc_time, Buffer,MaxEventDur)
        Nevent = size(Voc_transc_time,1);
        OutData.SpikeSUVoc = cell(Nevent,1);
        OutData.SpikeSUBSL = cell(Nevent,1);
        % reframe the extraction windows of each vocalization
        [EventOnset_time ,EventOffset_time] = reframe(Voc_transc_time, Buffer, MaxEventDur);
        % loading the single unit spike arrival times
        Spikes = load(InputFile, 'Spike_arrival_times');
        % loop through vocalizations and extract spike arrival times
        for vv=1:Nevent
            % Find the spike arrival times that are between the
            % requested times and center them to the onset of the
            % behavioral event
            OutData.SpikeSUVoc{vv} = Spikes.Spike_arrival_times(logical((Spikes.Spike_arrival_times>(EventOnset_time(vv)*10^3)) .* (Spikes.Spike_arrival_times<(EventOffset_time(vv)*10^3))))/10^3 - Voc_transc_time(vv,1);
            OutData.SpikeSUBSL{vv} = Spikes.Spike_arrival_times(logical((Spikes.Spike_arrival_times>(BSL_transc_time(vv,1)*10^3)) .* (Spikes.Spike_arrival_times<(BSL_transc_time(vv,2)*10^3))))/10^3 - Voc_transc_time(vv,1);
        end
    end



% Enlarge the window of extraction
    function [EventOnset_time, EventOffset_time] = reframe(OnsetOffset_time, Buffer, MaxEventDur)
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
end
