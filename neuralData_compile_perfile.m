function neuralData_compile_perfile(InputDataFile, OutputPath, NeuralBuffer)

% INPUT:
%       InputDataFile: path to a single unit mat file.

%       OutputPath: string that indicate where the data should be saved. If
%       not specified or left empty, the data will be saved in the folder
%       of the inputfile


% OUTPUT
% Struct with all the data


% Hard coded parameters for ploting of the spectrum in biosound
F_high_Raw = 50000;
F_high_Piezo = 10000;
Debug_Fig=0;

% Get the paths
[Path2Data, DataFile]=fileparts(InputDataFile);
PathParts = regexp(Path2Data, '/', 'split');
Loggers_dir = ['/' fullfile(PathParts{1:end-2})];
Audio_dir = ['/' fullfile(PathParts{1:end-4}, 'audio',PathParts{end-2})];
if nargin<2 || isempty(OutputPath)
    OutputPath = Loggers_dir;
end

% Get the date of the recording
Idx_ = strfind(DataFile, '_');
Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));

% on my MAC
% Loggers_dir = ['/' fullfile(PathParts{1:end-1},Date)];

% Get the tetrode ID
NeuralInputID{1} = DataFile(strfind(DataFile, 'TT')+2);
% Get the SS ID
NeuralInputID{2} = DataFile((Idx_(end)+1):end);
% Get the SS quality
NeuralInputID{3} = DataFile(strfind(DataFile, '_SS')+3);

% Get the subject ID
SubjectID = DataFile(1:5);

% Output
OutputDataFile = fullfile(OutputPath, sprintf('%s_%s_SS%s_%s-%s.mat', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2}));

% Loop through files to cut the vocalizations into analyzed snippets
% classify vocalizations according to type: Chirp or Trills
% fill in a who column
% fill in a what column (Chirp, Trill)
% fill in a duration

%% Initialize output varables
% get the number of files and fieldnames, to create cells from the
% following variables:
NExpe = 0;
% Initialize output matrix
DataDir = dir(fullfile(OutputPath, sprintf('%s_%s_*_SS%s_%s-%s.mat', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2})));
for ff=1:length(DataDir)
    Neuro = load(fullfile(DataDir(ff).folder, DataDir(ff).name));
    FieldNames = fieldnames(Neuro);
    for fi=1:length(FieldNames)
        if strcmp(FieldNames{fi}, 'Voc_NeuroSSU')
            NExpe = NExpe + 1;
        elseif strcmp(FieldNames{fi}, 'Behav_NeuroSSU')
            NExpe = NExpe + 1;
        end
    end
end
Duration = cell(1,NExpe); % Duration of the behavioral event in ms
DelayBefore = cell(1,NExpe); % Duration of no behavioral event before the onset in ms if no event before, Delay ms is considered as the safest estimate of the silence time before but here we are doing the assumption that nothing happened between call detected by voc_localize and voc_localize_operant
DelayAfter = cell(1,NExpe);% Duration of no behavioral event after the offset in ms if no event after, Delay ms is considered as the safest estimate of the silence time after but here we are doing the assumption that nothing happened between call detected by voc_localize and voc_localize_operant
VocOverlap = cell(1,NExpe); % Indicate if the vocalization overlaps with another one. 0: no overlap; 1: overlap and at least one overlaping vocalization is louder; 2: overlap and it is the loudest vocalization
VocWave = cell(1,NExpe);% Wave of the vocalization exactly extracted on Mic
VocPiezoWave = cell(1,NExpe);% Wave of the vocalization exactly extracted on Piezo
VocRank = cell(1,NExpe); % Rank of the call in the vocalization sequence
RewardTime = cell(1,NExpe); % Time of the reward alligned to voc onset
BioSound = cell(1,NExpe);%
BSLDuration = cell(1,NExpe);% Duration of the baseline sequence
SpikesArrivalTimes_Baseline = cell(1,NExpe); % Spike arrival time of the Baseline sequence
SpikesArrivalTimes_Behav = cell(1,NExpe);% Spike arrival time of the behavioral event
Who = cell(1,NExpe);% Identity of the performing bat (self or ID of the bat)
What = cell(1,NExpe);% Type of Behavior
ExpType = cell(1,NExpe);% Type of experiment (P=Play-Back, O=Operant conditioning, F=Free interactions)
NEvents = nan(1,NExpe);%number of behavioral event for each experiment
 % At the end concatenate them using [C{:}] and reshape([C{:}], NTot,1)

%% Load all behavior actions
NExpe = 0; % reinitialize the counter of experiments
DataDir = dir(fullfile(OutputPath, sprintf('%s_%s_*_SS%s_%s-%s.mat', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2})));
for ff=1:length(DataDir)
    Neuro = load(fullfile(DataDir(ff).folder, DataDir(ff).name));
    FieldNames = fieldnames(Neuro);
    Idx_2 = strfind(DataDir(ff).name,'_');
    ExpStartTime = DataDir(ff).name((Idx_2(2) +1):(Idx_2(3)-1));
    fprintf(1, 'Experiment time: %s\n',ExpStartTime) 
    % finding session type
    ParamFile = dir(fullfile(Audio_dir, sprintf('*%s_%s*param.txt', Date(3:end),ExpStartTime)));
    if contains(ParamFile.name, 'VocTrigger')
        SessionType_local = 'O';
    elseif contains(ParamFile.name, 'RecOnly')
        SessionType_local = 'F';
    end
    for fi=1:length(FieldNames)
        if strcmp(FieldNames{fi}, 'Voc_NeuroSSU')
            fprintf(1, '  Vocalization data\n')
            NExpe = NExpe + 1; % Increment the counter of experiments
            AudioDir_all = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData*.mat', Date(3:end), ExpStartTime))); % These are all the results of vocalization localization corresponding to that experiment time
            BeforeCuration = nan(length(AudioDir_all),1);
            for audiof=1:length(AudioDir_all)
                BeforeCuration(audiof)=(length(strfind(AudioDir_all(audiof).name, '_'))==2);
            end
            DataFileWho = AudioDir_all(logical(~BeforeCuration)); % These are the files obtained after manual curation
            Idx_ = strfind(DataFileWho(1).name,'_');
            Idxmat = strfind(DataFileWho(1).name,'.mat');
            Delay = str2double(DataFileWho(1).name((Idx_(end)+1):(Idxmat-1)));
            
            % find the logger number worn by the subject
            load(fullfile(DataFileWho(1).folder, DataFileWho(1).name), 'BatID','LoggerName');
            ALNum = contains(LoggerName, 'AL');
            SubjectLogs = cell2mat(BatID) == str2double(SubjectID);
            SubjectAL = LoggerName{find(ALNum .* SubjectLogs)}(3:end); %#ok<FNDSB>
            
            % find the vocalizations emitted by the vocalizer of interest
            DataFileBeforeCuration = AudioDir_all(logical(BeforeCuration)); % These are the files obtained before manual curation
            load(fullfile(DataFileBeforeCuration(1).folder, DataFileBeforeCuration(1).name),'Piezo_wave');
            Fns_AL = fieldnames(Piezo_wave);
            FocIndAudio = find(contains(Fns_AL, SubjectAL));
            
            
            % Loop through files to First pre-allocate space for output
            % variables
            NFiles = length(DataFileBeforeCuration);
            if NFiles~=length(DataFileWho)
                warning('Issue with the number of files before and after manual curation of the vocalizations (before/after who_calls')
                keyboard
            end
            
            VocCall = zeros(1, NFiles);
            Nseq = [0 nan(1, NFiles)];
            for nf = 1:NFiles
                if ~strcmp(DataFileWho(nf).name(27), '_') && ~strcmp(DataFileBeforeCuration(nf).name(27), '.') && (str2double(DataFileWho(nf).name(27))~=str2double(DataFileBeforeCuration(nf).name(27)))
                    warning('Issue with files ordering before and after manual curation of the vocalizations (before/after who_calls')
                    keyboard
                end
                load(fullfile(DataFileWho(nf).folder, DataFileWho(nf).name), 'IndVocStartRaw_merged');
                load(fullfile(DataFileBeforeCuration(nf).folder, DataFileBeforeCuration(nf).name), 'VocFilename','Old_vv_out_list');
                
                Nseq(nf+1) = length(VocFilename);
                if length(IndVocStartRaw_merged)>length(VocFilename) && exist('Old_vv_out_list', 'var') % This is a correction of the correction of merge patch in Who_Calls_Play_Less :-P
                    
                    load(fullfile(DataFileWho(nf).folder, DataFileWho(nf).name), 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged','IndVocStartRaw','IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo','IndVocStart_all', 'IndVocStop_all', 'RMSRatio_all', 'RMSDiff_all');
                    IndVocStartRaw_merged = IndVocStartRaw_merged(1:length(VocFilename));
                    IndVocStopRaw_merged = IndVocStopRaw_merged(1:length(VocFilename));
                    IndVocStartPiezo_merged = IndVocStartPiezo_merged(1:length(VocFilename));
                    IndVocStopPiezo_merged = IndVocStopPiezo_merged(1:length(VocFilename));
                    IndVocStartRaw = IndVocStartRaw(1:length(VocFilename));
                    IndVocStopRaw = IndVocStopRaw(1:length(VocFilename));
                    IndVocStartPiezo = IndVocStartPiezo(1:length(VocFilename));
                    IndVocStopPiezo = IndVocStopPiezo(1:length(VocFilename));
                    IndVocStart_all  = IndVocStart_all(1:length(VocFilename));
                    IndVocStop_all  = IndVocStop_all(1:length(VocFilename));
                    RMSRatio_all = RMSRatio_all(1:length(VocFilename));
                    RMSDiff_all= RMSDiff_all(1:length(VocFilename));
                    save(fullfile(DataFileWho(nf).folder, DataFileWho(nf).name), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged','IndVocStartRaw','IndVocStopRaw', 'IndVocStartPiezo', 'IndVocStopPiezo','IndVocStart_all', 'IndVocStop_all', 'RMSRatio_all', 'RMSDiff_all', '-append');
                    
                    clear IndVocStopRaw_merged IndVocStartPiezo_merged IndVocStopPiezo_merged IndVocStartRaw IndVocStopRaw IndVocStartPiezo IndVocStopPiezo IndVocStart_all IndVocStop_all RMSRatio_all RMSDiff_all
                end
             
                % Call sequences with identified vocalizations and re-order the
                % vocalization in chornological order (calculated in
                % cut_neural_data_voc)
                VocInd_Bin = logical([zeros(1,sum(Nseq(1:nf))) ~cellfun('isempty',IndVocStartRaw_merged)]);
                VocInd = Neuro.Voc_NeuroSSU.SortInd(VocInd_Bin) - sum(Nseq(1:nf));
                IndSort_local = find(~cellfun('isempty',IndVocStartRaw_merged));
                if any(~(VocInd == IndSort_local'))
                    warning('Looks like there was some re-ordering in time happening between sequences')
                end
                NV = length(VocInd);

                % Count the number of vocalization cuts for preallocation of space
                for vv=1:NV
                    for ll=1:length(IndVocStartRaw_merged{VocInd(vv)})
                        VocCall(nf) = VocCall(nf) + length(IndVocStartRaw_merged{VocInd(vv)}{ll});
                    end
                end
            end

            
            
            % Keep track of the number of behavioral event for each
            % experiment
            NEvents(NExpe) = sum(VocCall);
            
            % Now loop through files and calls and gather data
            Duration{NExpe} = nan(1,sum(VocCall)); % Duration of the vocalization in ms
            DelayBefore{NExpe} = nan(1,sum(VocCall)); % Duration of no behavioral event before the onset in ms
            DelayAfter{NExpe} = nan(1,sum(VocCall));% Duration of no behavioral event after the offset in ms
            VocOverlap{NExpe} = nan(1,sum(VocCall));% Signal the occurence of anoverlap with another vocalization from another bat
            VocWave{NExpe} = cell(1,sum(VocCall));% Wave of the vocalization exactly extracted on Mic
            VocPiezoWave{NExpe} = cell(1,sum(VocCall));% Wave of the vocalization exactly extracted on Piezo
            VocRank{NExpe} = cell(1,sum(VocCall));% Rank of the vocal element in the sequence of vocalization as first last or middle
            RewardTime{NExpe} = nan(1,sum(VocCall)); % time of the reward if call was rewarded
            BioSound{NExpe} = cell(2,sum(VocCall));% Biosound of the microphone on row1, biosound of the piezo on line 2
            BSLDuration{NExpe} = nan(1,sum(VocCall));% Duration of the baseline sequence
            SpikesArrivalTimes_Baseline{NExpe} = cell(1,sum(VocCall)); % Spike arrival time of the Baseline sequence
            SpikesArrivalTimes_Behav{NExpe} = cell(1,sum(VocCall));% Spike arrival time of the call event
            Who{NExpe} = cell(1,sum(VocCall));% Identity of the performing bat (self or ID of the bat)
            What{NExpe} = cell(1,sum(VocCall));% Type of Behavior
            ExpType{NExpe} = cell(1,sum(VocCall)); % Type of experiment (P=Play-Back, O=Operant conditioning, F=Free interactions)
            AllVocCall = VocCall;
            
            VocCall = zeros(1,NFiles);
            
            for nf = 1:NFiles
                if ~strcmp(DataFileWho(nf).name(27), '_') && ~strcmp(DataFileBeforeCuration(nf).name(27), '.') && (str2double(DataFileWho(nf).name(27))~=str2double(DataFileBeforeCuration(nf).name(27)))
                    warning('Issue with files ordering before and after manual curation of the vocalizations (before/after who_calls')
                    keyboard
                end
                load(fullfile(DataFileWho(nf).folder, DataFileWho(nf).name), 'IndVocStartRaw_merged','IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged','BioSoundCalls');
                load(fullfile(DataFileBeforeCuration(nf).folder, DataFileBeforeCuration(nf).name), 'FS','Piezo_wave','Raw_wave','Piezo_FS','Voc_transc_time_refined');
                % Call sequences with identified vocalizations and re-order the
                % vocalization in chornological order (calculated in
                % cut_neural_data_voc)
                VocInd_Bin = logical([zeros(1,sum(Nseq(1:nf))) ~cellfun('isempty',IndVocStartRaw_merged)]);
                VocInd = Neuro.Voc_NeuroSSU.SortInd(VocInd_Bin) - sum(Nseq(1:nf));
                IndSort_local = find(~cellfun('isempty',IndVocStartRaw_merged));
                if any(~(VocInd == IndSort_local'))
                    warning('Looks like there was some re-ordering in time happening between sequences')
                end
                NV = length(VocInd);
                Ncall = nan(NV,sum(ALNum));
                for vv=1:NV
                    for ll=1:length(IndVocStartRaw_merged{VocInd(vv)})
                        Ncall(vv,ll) = length(IndVocStartRaw_merged{VocInd(vv)}{ll});
                        if Ncall(vv,ll)
                            % Get the vector of all starts and stops of
                            % vocalizations detected in that sequence and in
                            % the previous or following sequences
                            if vv==NV
                                AllStarts = cell2mat(IndVocStartRaw_merged{VocInd(vv)}');
                            else
                                jj=0;
                                Next=[];
                                while (isempty(Next)) && (jj<(length(VocInd)-vv))
                                    jj=jj+1;
                                    % next onset times if any
                                    Next = cell2mat(IndVocStartRaw_merged{VocInd(vv+jj)}');
                                    % timedelay between this sequence begining and
                                    % the next sequence onset
                                    Delay2next = (Voc_transc_time_refined(VocInd(vv+jj),1) - Voc_transc_time_refined(VocInd(vv),1))*FS/1000; %Voc_transc_time_refined is in ms
                                    AllStarts = [cell2mat(IndVocStartRaw_merged{VocInd(vv)}') (Next+Delay2next)];
                                end
                                if isempty(Next)
                                    AllStarts = cell2mat(IndVocStartRaw_merged{VocInd(vv)}');
                                end
                            end
                            if vv>1
                                jj=0;
                                Previous = [];
                                while (isempty(Previous)) && ((jj+1)<vv)
                                    jj=jj+1;
                                    % previous offset times if any
                                    Previous = cell2mat(IndVocStopRaw_merged{VocInd(vv-jj)}');
                                    % timedelay between this sequence begining and
                                    % the previous sequence onset
                                    Delay2previous = (Voc_transc_time_refined(VocInd(vv),1) - Voc_transc_time_refined(VocInd(vv-jj),1))*FS/1000; %Voc_transc_time_refined is in ms
                                    AllStops = [cell2mat(IndVocStopRaw_merged{VocInd(vv)}') (Previous- Delay2previous)];
                                end
                                if isempty(Previous)
                                    AllStops = cell2mat(IndVocStopRaw_merged{VocInd(vv)}');
                                end
                            else
                                AllStops = cell2mat(IndVocStopRaw_merged{VocInd(vv)}');
                            end
                            for nn=1:Ncall(vv,ll)
                                VocCall(nf) = VocCall(nf)+1; % Increment the counter of vocalization events
                                % Save the vocalization Rank as first or last
                                if nn==1
                                    VocRank{NExpe}{sum(VocCall)} = 'first';
                                elseif nn==Ncall(vv,ll)
                                    VocRank{NExpe}{sum(VocCall)} = 'last';
                                else
                                    VocRank{NExpe}{sum(VocCall)} = 'middle';
                                end
                                % Save the type of experiment
                                ExpType{NExpe}{sum(VocCall)} = SessionType_local;
                                % Identify if the vocalizing bat was the one
                                % from which the neural data come from
                                if ll == FocIndAudio
                                    Who{NExpe}{sum(VocCall)} = 'self';
                                else
                                    AL_local = Fns_AL{ll};
                                    ALNum = contains(LoggerName, ['AL' AL_local(7:end)]);
                                    Who{NExpe}{sum(VocCall)} = num2str(BatID{ALNum});
                                end
                                % Get the duration of the vocalization in ms, the
                                % time preceding it since the last event and
                                % the time following until the next event
                                % %% INEED TO TAKE INTO ACCOUNT WHAT'S
                                % HAPPENING ON OTHER LOGGERS!
                                Duration{NExpe}(sum(VocCall)) = (IndVocStopRaw_merged{VocInd(vv)}{ll}(nn) -IndVocStartRaw_merged{VocInd(vv)}{ll}(nn))/FS*1000;
                                Durations2Preceding_events = IndVocStartRaw_merged{VocInd(vv)}{ll}(nn) - AllStops;
                                if sum(Durations2Preceding_events>0) % some calls happened just before
                                    DelayBefore{NExpe}(sum(VocCall)) = (min(Durations2Preceding_events(Durations2Preceding_events>0)))/FS*1000;
                                elseif VocInd(vv)==1 % This is the first vocalization of the recording session assume that nothing happened before
                                    DelayBefore{NExpe}(sum(VocCall)) = NeuralBuffer;
                                else % At best vocalizations before were for sure at a minimum of 200ms (merging time)
                                    DelayBefore{NExpe}(sum(VocCall)) = Delay;
                                end
                                Durations2Following_events = AllStarts - IndVocStopRaw_merged{VocInd(vv)}{ll}(nn);
                                if sum(Durations2Following_events>0)
                                    DelayAfter{NExpe}(sum(VocCall)) = (min(Durations2Following_events(Durations2Following_events>0)))/FS*1000;
                                elseif VocInd(vv) == max(VocInd) && nn == Ncall(vv,ll)
                                    DelayAfter{NExpe}(sum(VocCall)) = NeuralBuffer;
                                else
                                    DelayAfter{NExpe}(sum(VocCall)) = Delay;
                                end
                                % Any overlap with another event? 0: no
                                % overlap, 1: Overlap and this event is the
                                % quietest, 2: Overlap and this event is
                                % the loudest
                                Onset = IndVocStartRaw_merged{VocInd(vv)}{ll}(nn);
                                Offset = IndVocStopRaw_merged{VocInd(vv)}{ll}(nn);
                                Stops_in_vv = cell2mat(IndVocStopRaw_merged{VocInd(vv)}');
                                Starts_in_vv = cell2mat(IndVocStartRaw_merged{VocInd(vv)}');
                                Overlap_or_included = find((Onset>Starts_in_vv) .* (Onset<Stops_in_vv) + (Offset>Starts_in_vv) .* (Offset<Stops_in_vv));
                                Occlude = find((Onset<Starts_in_vv) .* (Offset>Stops_in_vv));
                                All_overlap_ind = [Overlap_or_included Occlude];
                                if ~isempty(All_overlap_ind)
                                    OnsetInd = find(Starts_in_vv == Onset);
                                    OffsetInd = find(Stops_in_vv == Offset);
                                    if OnsetInd~=OffsetInd
                                        warning('Unexpected error, these two haveto be the same')
                                        keyboard
                                    end
                                    RMS_overlap = nan(length(All_overlap_ind),1);
                                    for ovi=1:length(All_overlap_ind)
                                        VocCall_overlap = VocCall(nf) + (All_overlap_ind(ovi) - OnsetInd);
                                        RMS_overlap(ovi) = BioSoundCalls{VocCall_overlap,2}.rms;
                                    end
                                    if all(BioSoundCalls{VocCall(nf),2}.rms>RMS_overlap)
                                        VocOverlap{NExpe}(sum(VocCall)) = 2;
                                    else
                                        VocOverlap{NExpe}(sum(VocCall)) = 1;
                                    end
                                else
                                    VocOverlap{NExpe}(sum(VocCall)) = 0;
                                end
                                
                                % Duration of the baseline sequence
                                VocInd_Neuro = find(Neuro.Voc_NeuroSSU.SortInd == (VocInd(vv) + sum(Nseq(1:nf))));
                                BSLDuration{NExpe}(sum(VocCall)) = Neuro.Voc_NeuroSSU.BSL_transc_time_refined(VocInd_Neuro,2)-Neuro.Voc_NeuroSSU.BSL_transc_time_refined(VocInd_Neuro,1);

                                % Extract the sound of the microphone that
                                % correspond to the data with the same delay
                                % before/after as the neural response
                                IndOn = round(IndVocStartRaw_merged{VocInd(vv)}{ll}(nn) - NeuralBuffer*FS/1000);
                                IndOff = round(IndVocStopRaw_merged{VocInd(vv)}{ll}(nn) + NeuralBuffer*FS/1000);
                                DurWave_local = length(Raw_wave{VocInd(vv)});
                                WL = Raw_wave{VocInd(vv)}(max(1,IndOn):min(DurWave_local, IndOff));
                                if IndOn<1
                                    WL = [zeros(abs(IndOn),1); WL]; %#ok<AGROW>
                                end
                                if IndOff>DurWave_local
                                    WL = [WL ; zeros(IndOff-DurWave_local,1)]; %#ok<AGROW>
                                end
                                VocWave{NExpe}{sum(VocCall)} = WL; % contains the vocalization with the same delay before/after as the neural response

                                % Extract the sound of the audio-logger that
                                % correspond to the neural data (NeuralBuffer
                                % before/after vocalization onset/offset)
                                Piezo_samprate = Piezo_FS.(Fns_AL{ll})(VocInd(vv));
                                IndOn = round(IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn) - NeuralBuffer*Piezo_samprate/1000);
                                IndOff = round(IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn) + NeuralBuffer*Piezo_samprate/1000);
                                DurWave_local = length(Piezo_wave.(Fns_AL{ll}){VocInd(vv)});
                                WL = Piezo_wave.(Fns_AL{ll}){VocInd(vv)}(max(1,IndOn):min(DurWave_local, IndOff));
                                WL = reshape(WL,length(WL),1);
                                if IndOn<1
                                    WL = [zeros(abs(IndOn),1); WL]; %#ok<AGROW>
                                end
                                if IndOff>DurWave_local
                                    WL = [WL ; zeros(IndOff-DurWave_local,1)]; %#ok<AGROW>
                                end
                                VocPiezoWave{NExpe}{sum(VocCall)} = WL; % contains the vocalization with the same delay NeuralBuffer before/after as the neural response

                                % Get the biosound parameters
                                BioSound{NExpe}{1,sum(VocCall)} = BioSoundCalls{VocCall(nf),1};
                                BioSound{NExpe}{2,sum(VocCall)} = BioSoundCalls{VocCall(nf),2};

                                % Identify the type of call VocTr for a Trill,
                                % VocBa for a bark and VocUn for undefined Voc
                                What{NExpe}{sum(VocCall)} = ['Voc' BioSoundCalls{VocCall(nf),1}.type];
    %                             What{NExpe}{VocCall} = identify_CallType(Raw_wave{VocInd(vv)}(IndVocStartRaw_merged{VocInd(vv)}{ll}(nn):IndVocStopRaw_merged{VocInd(vv)}{ll}(nn)),Piezo_wave.(Fns_AL{ll}){VocInd(vv)}(IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn):IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn)));

                                % Finding the spikes that are during, before
                                % NeuroBuffer ms and after NeuroBuffer ms the vocalization 
                                IndSU01 = logical((Neuro.Voc_NeuroSSU.SpikeSUVoc{VocInd_Neuro}>(IndVocStartRaw_merged{VocInd(vv)}{ll}(nn)/FS*1000 - NeuralBuffer)) .* (Neuro.Voc_NeuroSSU.SpikeSUVoc{VocInd_Neuro}<(IndVocStopRaw_merged{VocInd(vv)}{ll}(nn)/FS*1000 + NeuralBuffer)));
                                % saving the spike arrival times in ms after
                                % centering them to the vocalization onset
                                SpikesArrivalTimes_Behav{NExpe}{sum(VocCall)} = (Neuro.Voc_NeuroSSU.SpikeSUVoc{VocInd_Neuro}(IndSU01) - IndVocStartRaw_merged{VocInd(vv)}{ll}(nn)/FS*1000)';
                                % Finding the spikes that correspond to the
                                % call sequence baseline
                                if ~isnan(BSLDuration{NExpe}(sum(VocCall))) % if isnan: No baseline period could be taken before the onset of the vocalization
                                    if nn==1 % all calls cut within the sequence have the same baseline, only calcuating if first call
                                        SpikesArrivalTimes_Baseline{NExpe}{sum(VocCall)} = Neuro.Voc_NeuroSSU.SpikeSUBSL{VocInd_Neuro};
                                    else
                                        SpikesArrivalTimes_Baseline{NExpe}{sum(VocCall)} = SpikesArrivalTimes_Baseline{NExpe}{sum(VocCall)-1};
                                    end
                                end
                                % Save the reward time in ms after
                                % centering them to the vocalization onset
                                RewardTime{NExpe}(sum(VocCall)) = Neuro.Voc_NeuroSSU.ReTime(VocInd_Neuro)- IndVocStartRaw_merged{VocInd(vv)}{ll}(nn)/FS*1000;

                                % Debug figure if requested
                                if Debug_Fig
                                    if sum(VocCall)==1
                                        FIG1 = figure();
                                    else
                                        clf(FIG1)
                                    end
    %                                 SpikesInd  = logical((SpikesArrivalTimes_Behav{NExpe}{VocCall}> -DelayBefore{NExpe}(VocCall)) .*  (SpikesArrivalTimes_Behav{NExpe}{VocCall}< (DelayAfter{NExpe}(VocCall))+ Duration{NExpe}(VocCall)));
                                    plotBiosoundAndSpikes(BioSound{NExpe}{2,sum(VocCall)}, F_high_Piezo, NeuralBuffer, NeuralBuffer,Duration{NExpe}(sum(VocCall)), SpikesArrivalTimes_Behav{NExpe}{sum(VocCall)},RewardTime{NExpe}(sum(VocCall)),0)
                                    suplabel(sprintf('Total Voc %d/%d File %d/%d Voc %d/%d', sum(VocCall),sum(AllVocCall),nf, NFiles,VocCall(nf), AllVocCall(nf)),'t');
                                    pause(1)
                                else
                                    fprintf(1,'     Total Voc %d/%d File %d/%d Voc %d/%d\n', sum(VocCall),sum(AllVocCall),nf, NFiles,VocCall(nf), AllVocCall(nf))
                                end
                            end


                        end
                    end
                end
            end
            
        elseif strcmp(FieldNames{fi}, 'Behav_NeuroSSU')
            fprintf(1, '  Other behavior data\n')
            NExpe = NExpe + 1; % Increment the counter of experiments
            
            % Keep track of the number of behavioral event for each file
            NBehav = length(Neuro.Behav_NeuroSSU.SpikeSUBehav);
            NEvents(NExpe) = NBehav;
            
            % Now loop through behaviors and gather data
            Duration{NExpe} = Neuro.Behav_NeuroSSU.Duration'; % Duration of the behavior in ms
            DelayBefore{NExpe} = nan(1,NBehav); % Duration of no behavioral event before the onset in ms
            DelayAfter{NExpe} = nan(1,NBehav);% Duration of no behavioral event after the offset in ms
            VocOverlap{NExpe} = nan(1,NBehav);% This stays empty
            SpikesArrivalTimes_Behav{NExpe} = Neuro.Behav_NeuroSSU.SpikeSUBehav';% Spike arrival time of the behavioral event
            Who{NExpe} = cell(1,NBehav);% Identity of the performing bat (self or ID of the bat)
            What{NExpe} = Neuro.Behav_What(Neuro.Behav_NeuroSSU.BehavIdxSSUAlive)';% Type of Behavior
            ExpType{NExpe} = cell(1,NBehav); % Type of experiment (P=Play-Back, O=Operant conditioning, F=Free interactions)
            VocWave{NExpe} = cell(1,NBehav);% This stays empty
            VocPiezoWave{NExpe} = cell(1,NBehav);% This stays empty
            VocRank{NExpe} = cell(1,NBehav);% This stays empty
            BioSound{NExpe} = cell(2,NBehav);% This stays empty
            BSLDuration{NExpe} = nan(1,NBehav);% This stays empty
            RewardTime{NExpe} = nan(1,NBehav);% This stays empty
            SpikesArrivalTimes_Baseline{NExpe} = cell(1,NBehav); % This stays empty
            
            for bb=1:NBehav
                % Identify if the vocalizing bat was the one
                % from which the neural data come from
                if Neuro.Behav_Who(Neuro.Behav_NeuroSSU.BehavIdxSSUAlive(bb)) == str2double(SubjectID)
                    Who{NExpe}{bb} = 'self';
                else
                    Who{NExpe}{bb} = num2str(Neuro.Behav_Who(Neuro.Behav_NeuroSSU.BehavIdxSSUAlive(bb)));
                end
                ExpType{NExpe}{bb} = 'F';
                VocRank{NExpe}{bb} = ' ';% This has to be filled with some strings (here space) to sort later, cannot just be left empty
            end 
        end
    end
end

% Save Data as a single structure
Duration = reshape([Duration{:}],1,sum(NEvents))'; % Duration of the behavioral event in ms
DelayBefore = reshape([DelayBefore{:}],1,sum(NEvents))'; % Duration of no behavioral event before the onset in ms
DelayAfter = reshape([DelayAfter{:}],1,sum(NEvents))'; % Duration of no behavioral event after the offset in ms
VocOverlap = reshape([VocOverlap{:}],1,sum(NEvents))';% Indicate if the vocalization overlaps with another one. 0: no overlap; 1: overlap and at least one overlaping vocalization is louder; 2: overlap and it is the loudest vocalization
VocWave = reshape([VocWave{:}],1,sum(NEvents))'; % Wave of the vocalization exactly extracted on Mic
VocPiezoWave = reshape([VocPiezoWave{:}],1,sum(NEvents))'; % Wave of the vocalization exactly extracted on Piezo
VocRank = reshape([VocRank{:}], 1,sum(NEvents))';% Rank of the vocal element in the sequence of vocalization
BioSound = reshape([BioSound{:}],2,sum(NEvents))'; %
BSLDuration = reshape([BSLDuration{:}],1,sum(NEvents))'; % Duration of the baseline sequence
RewardTime = reshape([RewardTime{:}],1,sum(NEvents))'; % Time of the reward for operant vocalizations
SpikesArrivalTimes_Baseline = reshape([SpikesArrivalTimes_Baseline{:}],1,sum(NEvents))';  % Spike arrival time of the Baseline sequence
SpikesArrivalTimes_Behav = reshape([SpikesArrivalTimes_Behav{:}],1,sum(NEvents))'; % Spike arrival time of the behavioral event
Who = reshape([Who{:}],1,sum(NEvents))'; % Identity of the performing bat (self or ID of the bat)
What = reshape([What{:}],1,sum(NEvents))'; % Type of Behavior
ExpType = reshape([ExpType{:}],1,sum(NEvents))'; 

if exist(OutputDataFile, 'file')
    save(OutputDataFile, 'Duration','DelayBefore','DelayAfter', 'VocOverlap', 'VocWave', 'VocPiezoWave', 'VocRank', 'BioSound','BSLDuration', 'SpikesArrivalTimes_Baseline','SpikesArrivalTimes_Behav','Who','What','ExpType','RewardTime','-append');
else
    save(OutputDataFile, 'Duration','DelayBefore','DelayAfter','VocOverlap', 'VocWave', 'VocPiezoWave', 'VocRank', 'BioSound','BSLDuration', 'SpikesArrivalTimes_Baseline','SpikesArrivalTimes_Behav','Who','What','ExpType','RewardTime');
end

%% Local function


function plotBiosoundAndSpikes(BiosoundObj, F_high, DelayBefore, DelayAfter,Duration, SpikeArrivalTimes, RewardTime,FormantPlot)
        if nargin<7
            FormantPlot=1;
        end
        % Plot the results of biosound calculations
        ss1=subplot(2,1,1);
        ColorCode = get(groot,'DefaultAxesColorOrder');
        DBNOISE =12;
        f_low = 0;
        logB = - 20*log10(abs(double(BiosoundObj.spectro)));
        maxB = max(max(logB));
        minB = maxB-DBNOISE;
        
        imagesc(double(BiosoundObj.to)*1000,double(BiosoundObj.fo),logB);          % to is in seconds
        axis xy;
        caxis('manual');
        caxis([minB maxB]);
        cmap = spec_cmap();
        colormap(cmap);
        %         colorbar()
        
        v_axis = axis;
        v_axis(3)=f_low;
        v_axis(4)=F_high;
        v_axis(1) = -DelayBefore;
        v_axis(2) = Duration+DelayAfter;
        axis(v_axis);
        xlabel('Time (ms)'), ylabel('Frequency');
        
        
        % Plot the fundamental and formants if they were calculated
        %     if double(BiosoundFi.sal)>MinSaliency
        Legend = {'F0' 'Formant1' 'Formant2' 'Formant3'};
        IndLegend = [];
        if ~isempty(double(BiosoundObj.f0))
            hold on
            plot(double(BiosoundObj.to)*1000,double(BiosoundObj.f0),'r-','LineWidth',2)
            IndLegend = [1 IndLegend];
        end
        if FormantPlot
            hold on
            plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F1),'Color',ColorCode(4,:),'LineWidth',2)
            hold on
            plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F2),'Color',ColorCode(2,:),'LineWidth',2)
            hold on
            if any(~isnan(double(BiosoundObj.F3)))
                plot(double(BiosoundObj.to)*1000,double(BiosoundObj.F3),'Color',ColorCode(3,:),'LineWidth',2)
                IndLegend = [IndLegend 2:4];
            else
                IndLegend = [IndLegend 2:3];
            end
        end
        legend(Legend(IndLegend), 'Location','SouthWest')
        
        legend('AutoUpdate', 'off')
        % plot the spikes
        for ss=1:length(SpikeArrivalTimes)
            hold on
            plot(SpikeArrivalTimes(ss).*ones(2,1), [F_high-500 F_high-1500], '-k', 'LineWidth',2)
            hold on
        end
        hold off
        
        % Plot the reward time if there is any
        hold on
        if ~isnan(RewardTime)
            plot(RewardTime*ones(2,1), [f_low F_high], 'g-', 'LineWidth',2)
        end
        hold off
        
        ss2=subplot(2,1,2);
        yyaxis left
        plot((1:length(double(BiosoundObj.sound)))/BiosoundObj.samprate*1000,double(BiosoundObj.sound), 'k-','LineWidth',2)
        hold on
        YLIM = get(ss2,'YLim');
        YLIM = max(abs(YLIM)).*[-1 1];
        set(ss2, 'YLim', YLIM)
        ylabel('Amplitude')
        SoundAmp = double(py.array.array('d', py.numpy.nditer(BiosoundObj.amp)));
        yyaxis right
        plot(double(BiosoundObj.tAmp)*1000,double(SoundAmp), 'r-', 'LineWidth',2)
        YLIM = get(ss2,'YLim');
        YLIM = max(abs(YLIM)).*[-1 1];
        set(ss2, 'YLim', YLIM)
        set(ss2, 'XLim', v_axis(1:2))
        xlabel('Time (ms)')
        ylabel('Enveloppe')
        title(sprintf('AmpPeriodicity = %.3f AmpPF = %.1f Hz',BiosoundObj.AmpPeriodP, BiosoundObj.AmpPeriodF))
        % Plot the reward time if there is any
        hold on
        if ~isnan(RewardTime)
%             legend('AutoUpdate', 'off')
            plot(RewardTime*ones(2,1), [-1 1], 'g-', 'LineWidth',2)
        end
        hold off
        
        
    end
end