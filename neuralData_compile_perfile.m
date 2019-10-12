function neuralData_compile_perfile(InputDataFile, OutputPath)

% INPUT:
%       InputDataFile: path to a single unit mat file.

%       OutputPath: string that indicate where the data should be saved. If
%       not specified or left empty, the data will be saved in the folder
%       of the inputfile


% OUTPUT
% Struct with all the data

% Get the paths
[Path2Data, DataFile]=fileparts(InputDataFile);
PathParts = regexp(Path2Data, '/', 'split');
% Loggers_dir = ['/' fullfile(PathParts{1:end-2})];
Loggers_dir = '/Users/elie/Documents/ManipBats/LMC/20190607';
if nargin<2 || isempty(OutputPath)
    OutputPath = Loggers_dir;
end

% Output
OutputDataFile = fullfile(Path2Data, [DataFile(1:end-3) '_FullDataSet.mat']);

% Get the date of the recording
Idx_ = strfind(DataFile, '_');
Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));

% Get the tetrode ID
NeuralInputID{1} = DataFile(strfind(DataFile, 'TT')+2);
% Get the SS ID
NeuralInputID{2} = DataFile((Idx_(end)+1):end);

% Get the subject ID
SubjectID = DataFile(1:5);

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
DataDir = dir(fullfile(OutputPath, sprintf('%s_%s_*_SSU%s-%s.mat', SubjectID, Date,NeuralInputID{1},NeuralInputID{2})));
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
Duration = cell(NExpe,1); % Duration of the behavioral event in ms
DelayBefore = cell(NExpe,1); % Duration of no behavioral event before the onset in ms
DelayAfter = cell(NExpe,1);% Duration of no behavioral event after the offset in ms
VocWave = cell(NExpe,1);% Wave of the vocalization exactly extracted on Mic
VocPiezoWave = cell(NExpe,1);% Wave of the vocalization exactly extracted on Piezo
BioSound = cell(NExpe,1);%
BSLDuration = cell(NExpe,1);% Duration of the baseline sequence
SpikesArrivalTimes_Baseline = cell(NExpe,1); % Spike arrival time of the Baseline sequence
SpikesArrivalTimes_Behav = cell(NExpe,1);% Spike arrival time of the behavioral event
Who = cell(NExpe,1);% Identity of the performing bat (self or ID of the bat)
What = cell(NExpe,1);% Type of Behavior
ExpType = cell(NExpe,1);
NEvents = nan(NExpe,1);
ExpStartTime_Old = Inf;
 % At the end concatenate them using [C{:}] and reshape([C{:}], NTot,1)

%% Load all behavior actions
NExpe = 0; % reinitialize the counter of experiments
DataDir = dir(fullfile(OutputPath, sprintf('%s_%s_*_SSU%s-%s.mat', SubjectID, Date,NeuralInputID{1},NeuralInputID{2})));
for ff=1:length(DataDir)
    Neuro = load(fullfile(DataDir(ff).folder, DataDir(ff).name));
    FieldNames = fieldnames(Neuro);
    for fi=1:length(FieldNames)
        if strcmp(FieldNames{fi}, 'Voc_NeuroSSU')
            NExpe = NExpe + 1; % Increment the counter of experiments
            Idx_2 = strfind(DataDir(ff).name,'_');
            ExpStartTime_new = DataDir(ff).name((Idx_2(2) +1):(Idx_2(3)-1));
            DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date(3:end), ExpStartTime_new)));
            Idx_ = strfind(DataFile.name,'_');
            Idxmat = strfind(DataFile.name,'.mat');
            Delay = str2double(DataFile.name((Idx_(end)+1):(Idxmat-1)));
            load(fullfile(DataFile.folder, DataFile.name), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'BatID','LoggerName','BioSoundCalls');
            load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date(3:end), ExpStartTime_new)), 'FS','Piezo_wave','Raw_wave','Piezo_FS');
            
            % find the logger number worn by the subject
            ALNum = contains(LoggerName, 'AL');
            SubjectLogs = cell2mat(BatID) == str2double(SubjectID);
            SubjectAL = LoggerName{find(ALNum .* SubjectLogs)}(3:end);
            
            % find the vocalizations emitted by the vocalizer of interest
            Fns_AL = fieldnames(Piezo_wave);
            FocIndAudio = find(contains(Fns_AL, SubjectAL));
            % Number of call sequences with identified vocalizations
            VocInd = find(~cellfun('isempty',IndVocStartRaw_merged));
            NV = length(VocInd);
            % Count the number of vocalization cuts for preallocation of space
            VocCall = 0;
            for vv=1:NV
                for ll=1:length(IndVocStartRaw_merged{VocInd(vv)})
                    VocCall = VocCall + length(IndVocStartRaw_merged{VocInd(vv)}{ll});
                end
            end
            % Keep track of the number of behavioral event for each file
            NEvents(NExpe) = VocCall;
            
            % Now loop through calls and gather data
            Duration{NExpe} = nan(VocCall,1); % Duration of the vocalization in ms
            DelayBefore{NExpe} = nan(VocCall,1); % Duration of no behavioral event before the onset in ms
            DelayAfter{NExpe} = nan(VocCall,1);% Duration of no behavioral event after the offset in ms
            VocWave{NExpe} = cell(VocCall,1);% Wave of the vocalization exactly extracted on Mic
            VocPiezoWave{NExpe} = cell(VocCall,1);% Wave of the vocalization exactly extracted on Piezo
            BioSound{NExpe} = cell(VocCall,2);
            BSLDuration{NExpe} = nan(1,VocCall);% Duration of the baseline sequence
            SpikesArrivalTimes_Baseline{NExpe} = cell(VocCall,1); % Spike arrival time of the Baseline sequence
            SpikesArrivalTimes_Behav{NExpe} = cell(VocCall,1);% Spike arrival time of the call event
            Who{NExpe} = cell(VocCall,1);% Identity of the performing bat (self or ID of the bat)
            What{NExpe} = cell(VocCall,1);% Type of Behavior
            ExpType{NExpe} = cell(VocCall,1); % Type of experiment (P=Play-Back, O=Operant conditioning, F=Free interactions)
            
            VocCall = 0;
            Ncall = nan(NV,1);
            for vv=1:NV
                for ll=1:length(IndVocStartRaw_merged{VocInd(vv)})
                    Ncall(vv) = length(IndVocStartRaw_merged{VocInd(vv)}{ll});
                    if Ncall(vv)
                        % Get the vector of all starts and stops of
                        % vocalizations detected in that sequence
                        AllStarts = cell2mat(IndVocStartRaw_merged{VocInd(vv)});
                        AllStops = cell2mat(IndVocStopRaw_merged{VocInd(vv)});
                        for nn=1:Ncall(vv)
                            VocCall = VocCall+1; % Increment the counter of vocalization events
                            
                            % Save the tye of experiment
                            if str2double(ExpStartTime_Old) < str2double(ExpStartTime_new) % I always do operant before the free session
                                ExpType{NExpe}{VocCall} = 'F';
                            else
                                ExpType{NExpe}{VocCall} = 'O';
                            end
                            % Identify if the vocalizing bat was the one
                            % from which the neural data come from
                            if ll == FocIndAudio
                                Who{NExpe}{VocCall} = 'self';
                            else
                                AL_local = Fns_AL(ll);
                                ALNum = contains(LoggerName, ['AL' AL_local(7:end)]);
                                Who{NExpe}{VocCall} = num2str(BatID{ALNum});
                            end
                            % Get the duration of the vocalization, the
                            % time preceding it since the last event and
                            % the time following until the next event
                            Duration{NExpe}(VocCall) = (IndVocStopRaw_merged{VocInd(vv)}{ll}(nn) -IndVocStartRaw_merged{VocInd(vv)}{ll}(nn))/FS*1000;
                            Durations2Preceding_events = IndVocStartRaw_merged{VocInd(vv)}{ll}(nn) - AllStops;
                            if sum(Durations2Preceding_events>0)
                                DelayBefore{NExpe}(VocCall) = (min(Durations2Preceding_events(Durations2Preceding_events>0)))/FS*1000;
                            else
                                DelayBefore{NExpe}(VocCall) = Delay;
                            end
                            Durations2Following_events = AllStarts - IndVocStopRaw_merged{VocInd(vv)}{ll}(nn);
                            if sum(Durations2Following_events>0)
                                DelayAfter{NExpe}(VocCall) = (min(Durations2Following_events(Durations2Following_events>0)))/FS*1000;
                            else
                                DelayAfter{NExpe}(VocCall) = Delay;
                            end
                            % Duration of the baseline sequence
%                             BSLDuration{NExpe}(VocCall) = BSL_transc_time_refined(VocInd(vv),2)-BSL_transc_time_refined(VocInd(vv),1);
                            
                            % Extract the sound of the microphone that
                            % correspond to the data
                            IndOn = IndVocStartRaw_merged{VocInd(vv)}{ll}(nn) - DelayBefore{NExpe}(VocCall)*FS/1000;
                            IndOff = IndVocStopRaw_merged{VocInd(vv)}{ll}(nn)+DelayAfter{NExpe}(VocCall)*FS/1000;
                            DurWave_local = length(Raw_wave{VocInd(vv)});
                            WL = Raw_wave{VocInd(vv)}(max(1,IndOn):min(DurWave_local, IndOff));
                            if IndOn<1
                                WL = [zeros(abs(IndOn),1); WL]; %#ok<AGROW>
                            end
                            if IndOff>DurWave_local
                                WL = [WL ; zeros(IndOff-DurWave_local,1)]; %#ok<AGROW>
                            end
                            VocWave{NExpe}{VocCall} = WL; % contains the vocalization with the same delay before/after as the neural response
                            
                            % Extract the sound of the audio-logger that
                            % correspond to the data
                            Piezo_samprate = Piezo_FS.(Fns_AL{ll})(VocInd(vv));
                            IndOn = round(IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn) - DelayBefore{NExpe}(VocCall)*Piezo_samprate/1000);
                            IndOff = round(IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn) + DelayAfter{NExpe}(VocCall)*Piezo_samprate/1000);
                            DurWave_local = length(Piezo_wave.(Fns_AL{ll}){VocInd(vv)});
                            WL = Piezo_wave.(Fns_AL{ll}){VocInd(vv)}(max(1,IndOn):min(DurWave_local, IndOff));
                            WL = reshape(WL,length(WL),1);
                            if IndOn<1
                                WL = [zeros(abs(IndOn),1); WL]; %#ok<AGROW>
                            end
                            if IndOff>DurWave_local
                                WL = [WL ; zeros(IndOff-DurWave_local,1)]; %#ok<AGROW>
                            end
                            VocPiezoWave{NExpe}{VocCall} = WL; % contains the vocalization with the same delay before/after as the neural response
                            
                            % Get the biosound parameters
                            BioSound{NExpe}{VocCall,1} = BioSoundCalls{VocCall,1};
                            BioSound{NExpe}{VocCall,2} = BioSoundCalls{VocCall,2};
                            
                            % Identify the type of call VocTr for a Trill,
                            % VocBa for a bark and VocUn for undefined Voc
                            What{NExpe}{VocCall} = ['Voc' BioSoundCalls{VocCall,1}.type];
%                             What{NExpe}{VocCall} = identify_CallType(Raw_wave{VocInd(vv)}(IndVocStartRaw_merged{VocInd(vv)}{ll}(nn):IndVocStopRaw_merged{VocInd(vv)}{ll}(nn)),Piezo_wave.(Fns_AL{ll}){VocInd(vv)}(IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn):IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn)));
                            
                            % Finding the spikes that are during, before
                            % and after the vocalization
                            IndSU01 = logical((Neuro.Voc_NeuroSSU.SpikeSUVoc{vv}>(IndVocStartRaw_merged{VocInd(vv)}{ll}(nn)/FS*1000 - DelayBefore{NExpe}(VocCall))) .* (Neuro.Voc_NeuroSSU.SpikeSUVoc{vv}<(IndVocStopRaw_merged{VocInd(vv)}{ll}(nn)/FS*1000 + DelayAfter{NExpe}(VocCall))));
                            % saving the spike arrival times in ms after
                            % centering them to the vocalization onset
                            SpikesArrivalTimes_Behav{NExpe}{VocCall} = (Neuro.Voc_NeuroSSU.SpikeSUVoc{VocInd(vv)}(IndSU01) - IndVocStartRaw_merged{VocInd(VocInd(vv))}{ll}(nn)/FS*1000)';
                            % Finding the spikes that correspond to the
                            % call sequence baseline
%                             if ~isnan(BSLDuration{NExpe}(VocCall)) % if isnan: No baseline period could be taken before the onset of the vocalization
%                                 if nn==1 % all calls cut within the sequence have the same baseline, only calcuating if first call
%                                     SpikesArrivalTimes_Baseline{NExpe}{VocCall} = Neuro.Voc_NeuroSSU.SpikeSUBSL{VocInd(vv)};
%                                 else
%                                     SpikesArrivalTimes_Baseline{NExpe}{VocCall} = SpikesArrivalTimes_Baseline{NExpe}{VocCall-1};
%                                 end
%                             end
                        end
                        
                        
                    end
                end
            end
            
        elseif strcmp(FieldNames{fi}, 'Behav_NeuroSSU')
            NExpe = NExpe + 1; % Increment the counter of experiments
            
            % Keep track of the number of behavioral event for each file
            NBehav = length(Neuro.Behav_NeuroSSU);
            NEvents(NExpe) = NBehav;
            
            % Now loop through behaviors and gather data
            Duration{NExpe} = Neuro.Behav_NeuroSSU.Duration; % Duration of the behavior in ms
            DelayBefore{NExpe} = nan(NBehav,1); % Duration of no behavioral event before the onset in ms
            DelayAfter{NExpe} = nan(NBehav,1);% Duration of no behavioral event after the offset in ms
            SpikesArrivalTimes_Behav{NExpe} = Neuro.Behav_NeuroSSU.SpikeSUBehav;% Spike arrival time of the behavioral event
            Who{NExpe} = cell(NBehav,1);% Identity of the performing bat (self or ID of the bat)
            What{NExpe} = Neuro.Behav_What(Neuro.Behav_NeuroSSU.BehavIdxSSUAlive);% Type of Behavior
            ExpType{NExpe} = cell(NBehav,1); % Type of experiment (P=Play-Back, O=Operant conditioning, F=Free interactions)
            ExpType{NExpe}{:} = deal('F');
            
            for bb=1:NBehav
                % Identify if the vocalizing bat was the one
                % from which the neural data come from
                if Neuro.Behav_Who{Neuro.Behav_NeuroSSU.BehavIdxSSUAlive(bb)} == str2double(SubjectID)
                    Who{NExpe}{bb} = 'self';
                else
                    Who{NExpe}{bb} = num2str(Neuro.Behav_Who);
                end
            end 
        end
        ExpStartTime_Old = ExpStartTime_new;
    end
end

% Save Data as a single structure
Duration = reshape([Duration{:}],sum(NEvents),1); % Duration of the behavioral event in ms
DelayBefore = reshape([DelayBefore{:}],sum(NEvents),1); % Duration of no behavioral event before the onset in ms
DelayAfter = reshape([DelayAfter{:}],sum(NEvents),1); % Duration of no behavioral event after the offset in ms
VocWave = reshape([VocWave{:}],sum(NEvents),1); % Wave of the vocalization exactly extracted on Mic
VocPiezoWave = reshape([VocPiezoWave{:}],sum(NEvents),1); % Wave of the vocalization exactly extracted on Piezo
BioSound = reshape([BioSound{:}],sum(NEvents),2); %
BSLDuration = reshape([BSLDuration{:}],sum(NEvents),1); % Duration of the baseline sequence
SpikesArrivalTimes_Baseline = reshape([SpikesArrivalTimes_Baseline{:}],sum(NEvents),1);  % Spike arrival time of the Baseline sequence
SpikesArrivalTimes_Behav = reshape([SpikesArrivalTimes_Behav{:}],sum(NEvents),1); % Spike arrival time of the behavioral event
Who = reshape([Who{:}],sum(NEvents),1); % Identity of the performing bat (self or ID of the bat)
What = reshape([What{:}],sum(NEvents),1); % Type of Behavior
ExpType = reshape([ExpType{:}],sum(NEvents),1); 

save(OutputDataFile, 'Duration','DelayBefore','DelayAfter', 'VocWave', 'VocPiezoWave', 'BioSound','BSLDuration', 'SpikesArrivalTimes_Baseline','SpikesArrivalTimes_Behav','Who','What','ExpType')

end