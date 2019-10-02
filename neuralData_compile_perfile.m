function neuralData_compile_perfile(InputDataFile, OutputPath, MaxDur)

% INPUT:
%       InputDataFile: path to a single unit mat file.

%       OutputPath: string that indicate where the data should be saved. If
%       not specified or left empty, the data will be saved in the folder
%       of the inputfile
if nargin<3
    MaxDur = 500; % Duration in ms for estimating the rate
end

% Get the paths
[Path2Data, DataFile]=fileparts(InputDataFile);
PathParts = regexp(Path2Data, '/', 'split');
Loggers_dir = ['/' fullfile(PathParts{1:end-2})];
if nargin<2 || isempty(OutputPath)
    OutputPath = Loggers_dir;
end

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

%% Load the other behavior actions
% Initialize output matrix
DataDir = dir(fullfile(OutputPath, sprintf('%s_%s_*_SSU%s-%s.mat', SubjectID, Date,NeuralInputID{1},NeuralInputID{2})));
for ff=1:length(DataDir)
    Neuro = load(fullfile(DataDir(ff).folder, DataDir(ff).name));
    FieldNames = fieldnames(Neuro);
    for fi=1:length(FieldNames)
        if strcmp(FieldNames{fi}, 'Voc_NeuroSSU')
            Idx_2 = strfind(DataDir(ff).name,'_');
            ExpStartTime = DataDir(ff).name((Idx_2(2) +1):(Idx_2(3)-1));
            DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date(3:end), ExpStartTime)));
            load(fullfile(DataFile.folder, DataFile.name), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged', 'BatID','LoggerName');
            load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date(3:end), ExpStartTime)), 'FS','Piezo_wave','Raw_wave');
            
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
            % Now loop through calls and gather data
            Duration = nan(VocCall,1); % Duration of the vocalization in ms
            DelayBefore = nan(VocCall,1); % Duration of no behavioral event before the onset in ms
            DelayAfter = nan(VocCall,1);% Duration of no behavioral event after the offset in ms
            VocWave = cell(VocCall,1);% Wave of the vocalization exactly extracted on Mic
            VocPiezoWave = cell(VocCall,1);% Wave of the vocalization exactly extracted on Piezo
            BSLDuration = nan(1,VocCall);% Duration of the baseline sequence
            SpikesArrivalTimes_Baseline = cell(VocCall,1); % Spike arrival time of the Baseline sequence
            SpikesArrivalTimes_Call = cell(VocCall,1);% Spike arrival time of the call event
            Who = cell(VocCall,1);% Identity of the performing bat (self or ID of the bat)
            What = cell(VocCall,1);% Type of Behavior
            
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
                            
                            % Identify if the vocalizing bat was the one
                            % from which the neural data come from
                            if ll == FocIndAudio
                                Who{VocCall} = 'self';
                            else
                                AL_local = Fns_AL(ll);
                                ALNum = contains(LoggerName, ['AL' AL_local(7:end)]);
                                Who{VocCall} = num2str(BatID{ALNum});
                            end
                            % Get the duration of the vocalization, the
                            % time preceding it since the last event and
                            % the time following until the next event
                            Duration(VocCall) = (IndVocStopRaw_merged{VocInd(vv)}{ll}(nn) -IndVocStartRaw_merged{VocInd(vv)}{ll}(nn))/FS*1000;
                            Durations2Preceding_events = IndVocStartRaw_merged{VocInd(vv)}{ll}(nn) - AllStops;
                            DelayBefore(VocCall) = min(Durations2Preceding_events(Durations2Preceding_events>0));
                            if isempty(DelayBefore(VocCall))
                                DelayBefore(VocCall) = Delay;
                            end
                            Durations2Following_events = AllStarts - IndVocStopRaw_merged{VocInd(vv)}{ll}(nn);
                            DelayAfter(VocCall) = min(Durations2Following_events(Durations2Following_events>0));
                            if isempty(DelayAfter(VocCall))
                                DelayAfter(VocCall) = Delay;
                            end
                            % Duration of the baseline sequence
                            BSLDuration(VocCall) = BSL_transc_time_refined(VocInd(vv),2)-BSL_transc_time_refined(VocInd(vv),1);
                            
                            % Extract the sound of the microphone that
                            % correspond to the data
                            IndOn = IndVocStartRaw_merged{VocInd(vv)}{ll}(nn) - DelayBefore(VocCall)*FS/1000;
                            IndOff = IndVocStopRaw_merged{VocInd(vv)}{ll}(nn)+DelayAfter(VocCall)*FS/1000;
                            DurWave_local = length(Raw_wave{VocInd(vv)});
                            WL = Raw_wave{VocInd(vv)}(max(1,IndOn):min(DurWave_local, IndOff));
                            if IndOn<1
                                WL = [zeros(abs(IndOn),1); WL]; %#ok<AGROW>
                            end
                            if IndOff>DurWave_local
                                WL = [WL ; zeros(IndOff-DurWave_local,1)]; %#ok<AGROW>
                            end
                            VocWave{VocCall} = WL; % contains the vocalization with the same delay before/after as the neural response
                            
                            % Extract the sound of the audio-logger that
                            % correspond to the data
                            IndOn = IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn) - DelayBefore(VocCall)*FS_Piezo/1000;
                            IndOff = IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn) + DelayAfter(VocCall)*FS_Piezo/1000;
                            DurWave_local = length(Piezo_wave.(Fns_AL{ll}){VocInd(vv)});
                            WL = Piezo_wave.(Fns_AL{ll}){VocInd(vv)}(max(1,IndOn):min(DurWave_local, IndOff));
                            if IndOn<1
                                WL = [zeros(abs(IndOn),1); WL]; %#ok<AGROW>
                            end
                            if IndOff>DurWave_local
                                WL = [WL ; zeros(IndOff-DurWave_local,1)]; %#ok<AGROW>
                            end
                            VocPiezoWave{VocCall} = WL; % contains the vocalization with the same delay before/after as the neural response
                            
                            % Identify the type of call VocT for a Trill,
                            % VocC for a chirp and VocU for undefined Voc
                            What{VocCall} = identify_CallType(Raw_wave{VocInd(vv)}(IndVocStartRaw_merged{VocInd(vv)}{ll}(nn):IndVocStopRaw_merged{VocInd(vv)}{ll}(nn)),Piezo_wave.(Fns_AL{ll}){VocInd(vv)}(IndVocStartPiezo_merged{VocInd(vv)}{ll}(nn):IndVocStopPiezo_merged{VocInd(vv)}{ll}(nn)));
                            
                            % Finding the spikes that are during, before
                            % and after the vocalization
                            IndSU01 = logical((Neuro.Voc_NeuroSSU.SpikeSUVoc{vv}>(IndVocStartRaw_merged{VocInd(vv)}{ll}(nn)/FS*1000 - DelayBefore(VocCall))) .* (Neuro.Voc_NeuroSSU.SpikeSUVoc{vv}<(IndVocStopRaw_merged{VocInd(vv)}{ll}(nn)/FS*1000 + DelayAfter(VocCall))));
                            % saving the spike arrival times in ms after
                            % centering them to the vocalization onset
                            SpikesArrivalTimes_Call{VocCall} = (Neuro.Voc_NeuroSSU.SpikeSUVoc{VocInd(vv)}(IndSU01) - IndVocStartRaw_merged{VocInd(VocInd(vv))}{ll}(nn)/FS*1000)';
                            % Finding the spikes that correspond to the
                            % call sequence baseline
                            if ~isnan(BSLDuration(VocCall)) % if isnan: No baseline period could be taken before the onset of the vocalization
                                if nn==1 % all calls cut within the sequence have the same baseline, only calcuating if first call
                                    SpikesArrivalTimes_Baseline{VocCall} = Neuro.Voc_NeuroSSU.SpikeSUBSL{VocInd(vv)};
                                else
                                    SpikesArrivalTimes_Baseline{VocCall} = SpikesArrivalTimes_Baseline{VocCall-1};
                                end
                            end
                        end
                        
                        
                    end
                end
            end
            
            
        end
    end
end
Num_Slots(bb) = sum(floor(Dur_local/MaxDur)) + sum(mod(Dur_local, MaxDur)>0);

%% INTERNAL FUNCTION
function [What] = identify_CallType(MicWave, PiezoWave)
    

end
end