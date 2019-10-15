function get_logger_data_behav(Audio_dir, Loggers_dir, Date, ExpStartTime)
%% This function extract the onset/offset times of behaviors annotated during RecOnly using Matlab
% First it loads the files generated by code_social_behav.m (operant_bats
% library).
% Second it converts all the matlab computer time stamps of behavioral
% events to transceiver time using the output of
% align_soundmexAudio_2_logger applied on the RecOnly session.
% Third, it organizes data such as to return a cell array of onset/offset
% times of behaviors (AllActions_Time) where each cell corresponds to one type of behavior
% referenced in UActionText and where the first column of the N*2 matrix
% contained in that cell corresponds to onset of the behavior and second column to offset,
% each row being one instance of the behavior. The offset is calculated as 
% time when the key was hit - 1s to account for the delay to manually
% detect.
%  If the behavior is ponctual, then it returns only one
% time value (one column vector per cell). The code also returns the ID of
% the bat performing the behavior (AllActions_ID) for these behaviors where the iD can
% easily be retrieved (chewing, licking, quiet, teeth cleaning)
% Data are saved in an output file in the loggers directory, the name of
% the file following this rule:
% [ExpDate]_[ExpStartTime]_BehavExtractData.mat


%% Load the manual coding of behavior
InputFiles = dir(fullfile(Audio_dir,'*RecOnly_behav.txt'));
NF = length(InputFiles);
Type = cell(NF,1); % Type of annotation (vocalization, licking_start, chewing....)
Stamp = cell(NF,1); % Time stamp in computer time of the form yyyymmddThhmmssmmm
BatID = cell(NF,1); % ID of the bat for licking, quiet, chewing

% Go through files and gather data
for ff = 1:NF
    DataFile = fopen(fullfile(InputFiles(ff).folder, InputFiles(ff).name));
    textscan(DataFile, '%s\t%s\t%s\n',1);
    Data = textscan(DataFile, '%s\t%s\t%f');
    fclose(DataFile);
    LastRow = find(strcmp(Data{:,1},'stop'));
    if isempty(LastRow)
        LastRow = length(Data{1})+1;
    end
    Type{ff} = Data{1}(1:(LastRow-1));
    Stamp{ff} = Data{2}(1:(LastRow-1));
    IndRec = strfind(InputFiles(ff).name, 'RecOnly');
    Batname = InputFiles(ff).name((IndRec-3):(IndRec-2));
    if strcmp(Batname, 'Ha')
        BatID{ff} = 59882*ones(LastRow-1,1);
    elseif strcmp(Batname, 'Ho')
        BatID{ff} = 11689*ones(LastRow-1,1);
    elseif strcmp(Batname, 'Co')
        BatID{ff} = 59834*ones(LastRow-1,1);
    elseif strcmp(Batname, 'Ed')
        BatID{ff} = 65701*ones(LastRow-1,1);
    else
        BatIDLocal = input(sprintf('Imposible to identify Bat ID, please enter Chip ID for file %s', InputFiles(ff).name));
        BatID{ff} = BatIDLocal*ones(flip(size(1:(LastRow-1))));
    end
end
Stamp = vertcat(Stamp{:});
BatID = vertcat(BatID{:});
Type = vertcat(Type{:});

%% Convert computer time stamps to transceiver time
% Covert computer time stamps to sample stamps on the sound card
FS = 192000; % Sample frequency of the sound card (Not nice to set that hard here, but better than taking the time to load a wavfile to get the FS)
FormatIn = 'yyyymmddTHHMMSSFFF';% Format of the time in the raw files
% get the onset of the task
ParamFile = dir(fullfile(Audio_dir, sprintf('*%s_%s*RecOnly_param.txt', Date, ExpStartTime)));
if length(ParamFile)>1
    fprintf(1,'Several RecOnly sessions occured, which one do you want to look at?\n')
    for ee=1:length(ParamFile)
        fprintf(1,'%d: %s\n',ee,ParamFile(ee).name);
    end
    E_input = input('Your Choice:\n');
    ParamFile = ParamFile(E_input);
end
ParamF = fopen(fullfile(ParamFile.folder, ParamFile.name));
ParamData = textscan(ParamF, '%s', 'Delimiter','\n');
fclose(ParamF);
IndLine = find(contains(ParamData{1}, 'Task starts at '));
IndChar = length('Task starts at ')+1;
FirstFile = ParamData{1}{IndLine}(IndChar:end); %#ok<FNDSB>

% get the file changes time
EventFile = dir(fullfile(Audio_dir, sprintf('*%s_%s*RecOnly_events.txt', Date, ExpStartTime)));
if length(EventFile)>1
    fprintf(1,'Several RecOnly sessions occured, which one do you want to look at?\n')
    for ee=1:length(EventFile)
        fprintf(1,'%d: %s\n',ee,EventFile(ee).name);
    end
    E_input = input('Your Choice:\n');
    EventFile = EventFile(E_input);
end
EventF = fopen(fullfile(EventFile.folder, EventFile.name));
textscan(EventF,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',1);
EventData = textscan(EventF, '%s\t%f\t%s\t%s\t%d\t%d\t%f');
fclose(EventF);
ChangeFileInd = find(contains(EventData{4}, 'ChangeFile'));
ChangeFileStamp = EventData{1}(ChangeFileInd); %#ok<FNDSB> % This is the time stamp in computer time when a new microphone file was started

% Load the pulse times and samples
TTL_dir = dir(fullfile(Audio_dir,sprintf('%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
if isempty(TTL_dir)
    fprintf(1,'*** Alligning TTL pulses for the RecOnly session ***\n');
    align_soundmexAudio_2_logger(Audio_dir, Loggers_dir, ExpStartTime,'TTL_pulse_generator','Avisoft','Method','risefall', 'Session_strings', {'rec only start', 'rec only stop'});
    TTL_dir = dir(fullfile(Audio_dir,sprintf('%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
else
    fprintf(1,'*** ALREADY DONE: Alligning TTL pulses for the RecOnly session ***\n');
end
TTL = load(fullfile(TTL_dir.folder, TTL_dir.name));

 % loop through behavioral events and get the time in audio samples
 NE = length(Stamp);
 Event_AudioFileID = nan(NE,1);% File # during which each behavioral event happened
 Event_AudioStamp = nan(NE,1);% Time in sample of audiofile of each behavioral event
 Event_TranscTime = nan(NE,1); % Time in ms of each behavioral event
 for ee=1:NE
     LocalStamp = datevec(Stamp{ee}, FormatIn);
     for fc=1:length(ChangeFileStamp)
         FSt = datevec(ChangeFileStamp{fc}, FormatIn);
         Elapsed=etime(LocalStamp, FSt);
         if Elapsed<0
             Event_AudioFileID(ee) = fc;
             if fc==1 % This event happened during the first file
                 FSt = datevec(FirstFile, FormatIn);
             else
                FSt = datevec(ChangeFileStamp{fc-1}, FormatIn);
             end
             Elapsed=etime(LocalStamp, FSt);
             Event_AudioStamp(ee) = round(Elapsed*FS);
             break
         end
     end
     if (fc==length(ChangeFileStamp)) && isnan(Event_AudioFileID(ee))
         Event_AudioFileID(ee) = fc+1;
         FSt = datevec(ChangeFileStamp{fc}, FormatIn);
         Elapsed=etime(LocalStamp, FSt);
         Event_AudioStamp(ee) = round(Elapsed*FS);
     end
     
     % Extract the transceiver time
     % zscore the sample stamps
     TTL_idx = find(unique(TTL.File_number) == Event_AudioFileID(ee));
     Event_AudioStamp_zs = (Event_AudioStamp(ee) - TTL.Mean_std_Pulse_samp_audio(TTL_idx,1))/TTL.Mean_std_Pulse_samp_audio(TTL_idx,2);
     % calculate the transceiver times
     Event_TranscTime(ee) = TTL.Mean_std_Pulse_TimeStamp_Transc(TTL_idx,2) .* polyval(TTL.Slope_and_intercept{TTL_idx},Event_AudioStamp_zs,[], TTL.Mean_std_x{TTL_idx}) + TTL.Mean_std_Pulse_TimeStamp_Transc(TTL_idx,1);
 end


%% Identify individual actions and gather onset and offset times of each behavior
% Individualized actions:
IndivActions = {'chewing' 'food-in-mouth' 'licking' 'teeth-cleaning' 'quiet'};
UType_all = unique(Type); % unique labels
% Isolate actions items
StartInd= strfind(UType_all, 'start');
StopInd= strfind(UType_all, 'stop');
UActionText = cell(size(UType_all));
for ii=1:length(UType_all)
    if ~isempty(StartInd{ii})
        UActionText{ii} = UType_all{ii}(1:(StartInd{ii}-2));
    elseif ~isempty(StopInd{ii})
        UActionText{ii} = UType_all{ii}(1:(StopInd{ii}-2));
    else
        UActionText{ii} = UType_all{ii};
    end
end
UActionText = unique(UActionText);% different types of behaviors

% Extract the behavior onset/offset times for each bat
% % Find the reference time to plot the ethogram
% RefTime = Eventfile.event_timestamps_usec(find(contains(Eventfile.event_types_and_details, 'Started recording'),1))*10^-3;
NAction = length(UActionText);
AllActions_Time = cell(NAction,1);% onset/offset Time in ms of each long behavioral event or just time in msec for ponctual events at the scale of annotations (vocalizations for instance)
AllActions_ID = cell(NAction,1);

% Get the indices of all the starting points
IndStart = find(contains(Type, 'start'));
% Get the indices of all the stoping points
IndStop = find(contains(Type, 'stop'));
% Fig=figure();
% ColorCode = get(groot,'DefaultAxesColorOrder');
% ColorCode = [ColorCode(1:6,:) ; 0 0 1; ColorCode(7,:); 0.85 0.6940 0.556; 0 0 0; 1 0 0; 0.301 0.078 0.741; ];
% loop through each action
for aa=1:NAction
    fprintf('extracting data of %s\n',UActionText{aa})
    IndAction = find(contains(Type, UActionText{aa}));
    IndActionStart = intersect(IndAction, IndStart);
    if isempty(IndActionStart) % This action is ponctual
        fprintf('Saving ponctual events\n')
        AllActions_Time{aa} = Event_TranscTime(IndAction);
        % save the bat ID for that behavior if it is individualized
        if sum(strcmp(IndivActions, UActionText{aa}))
            AllActions_ID{aa} = BatID(IndAction);
        else
            AllActions_ID{aa} = nan(size(IndAction));
        end
%         hold on
%         plot((AllActions_Time{aa}-RefTime)*10^-3, aa, '*','Color', ColorCode(aa,:), 'MarkerSize', 5)
    else
        fprintf('Finding start and end of each event\n')
        IndActionStop = intersect(IndAction, IndStop);
        AllActions_Time{aa} = nan(max(length(IndActionStop), length(IndActionStart)),2);
        % allocate space for the bat ID for that behavior
        AllActions_ID{aa} = nan(max(length(IndActionStop), length(IndActionStart)),1);
        
        % loop through action starts and figure out if a stop is closer to
        % it than the next start
        e_count = 0;
        for ii=1:length(IndActionStart)
            fprintf('potential start %d/%d\n',ii,length(IndActionStart));
            % save the bat ID for that behavior
            BatID_local = BatID(IndActionStart(ii));
            % get the stop closer to the start from the same bat
            IndActionStopBat = intersect(IndActionStop, find(BatID == BatID_local));
            mm = find(min(IndActionStopBat - IndActionStart(ii)));
            if (IndActionStopBat(mm) - IndActionStart(ii)) <0
                error('Error in parsing the action, the stop happens before the start\n')
            end
            if isempty(IndActionStopBat)  % This start misses a stop, consider the next action as
                % being a potential stop
                IndActionStop_local = IndActionStart(ii)+1;
                if diff(Event_TranscTime(IndActionStart(ii):IndActionStop_local))>1000% only keep this behavioral event if the next behavior happens 1 sec later
                    e_count= e_count+1;
                    AllActions_Time{aa}(e_count,1) = Event_TranscTime(IndActionStart(ii));
                    AllActions_Time{aa}(e_count,2) = Event_TranscTime(IndActionStop_local) - 1000;
                    % save the bat ID for that behavior if it is individualized
                    if sum(strcmp(IndivActions, UActionText{aa}))
                        AllActions_ID{aa}(e_count) = BatID_local;
                    end
                end
            elseif ~isempty(IndActionStopBat) && ((ii==length(IndActionStart)) || (IndActionStopBat(mm) - IndActionStart(ii)) < diff(IndActionStart(ii:(ii+1))))  % this stop happens before another start, all is good save that pair of event
                e_count= e_count+1;
                AllActions_Time{aa}(e_count,1) = Event_TranscTime(IndActionStart(ii));
                AllActions_Time{aa}(e_count,2) = Event_TranscTime(IndActionStopBat(mm));
                % save the bat ID for that behavior if it is individualized
                if sum(strcmp(IndivActions, UActionText{aa}))
                    AllActions_ID{aa}(e_count) = BatID_local;
                end
                if ii<length(IndActionStart) && length(IndActionStop)>1
                    % Find the next start
                    nn = find(min(IndActionStart(ii+1:end) - IndActionStopBat(mm)));
                    kk = find(IndActionStop == IndActionStopBat(mm)); %Position of that particular stop in the list of potential stops from all bats
                    if (IndActionStart(ii+nn) - IndActionStopBat(mm)) < diff(IndActionStop(kk:(kk+1)))
                        % The next start happens before another stop, all is
                        % good, proceed to the next start and supress that stop
                        IndActionStop(kk) = [];
                    else % There are several stops in a row, get rid of them up to the next start
                        SupInd = find(diff(IndActionStop(kk:end))<(IndActionStart(ii+nn) - IndActionStop(kk)));
                        IndActionStop(kk : (kk+ SupInd)) = [];
                    end
                else %no more stops available stop here
                    break
                end
            end
%           % plot the action
%           hold on
%           plot((AllActions_Time{aa}(e_count,:)-RefTime)*10^-3, aa*ones(2,1), '-', 'Color', ColorCode(aa,:), 'LineWidth',2);
        end
        AllActions_Time{aa} = AllActions_Time{aa}(1:e_count,:);
        AllActions_ID{aa} = AllActions_ID{aa}(1:e_count);
    end
end

% Save data
save(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData.mat', Date, ExpStartTime)), 'AllActions_Time', 'AllActions_ID',  'UActionText');

end
    


