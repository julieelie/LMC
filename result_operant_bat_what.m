function result_operant_bat_what(Path2ParamFile, Path2RecordingTable, Logger_dir)
addpath(genpath('/Users/elie/Documents/CODE/LMC'))
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'))
addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))
TranscExtract = 1; % set to 1 to extract logger data and transceiver time
ForceExtract = 0; % set to 1 to redo the extraction of loggers otherwise the calculations will use the previous extraction data
ForceAllign = 0; % In case the TTL pulses allignment was already done but you want to do it again, set to 1
ForceVocExt3 = 0; % Used to patch the error in voc_localize_operant with call already detected issue
ForceVocExt1 = 0; % In case the extraction of vocalizations that triggered rewarding system was already done but you want to do it again set to 1
ForceVocExt2 = 0; % In case the extraction of vocalizations that triggered rewarding system was already done but you want to do it again set to 1
ForceWhoID = 0; % In case the identification of bats was already done but you want to re-do it again
ForceWhat = 1; % In case running biosound was already done but you want to re-do it
PlotIndivFile = 0; % Set to 1 to plot the sound pressure waveforms of individual detected vocalizations
close all
% Get the recording data
[AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
Date = DataFile(6:11);
WavFileStruc = dir(fullfile(AudioDataPath, [DataFile(1:16) '*mic*.wav']));

% Get the sound snippets from the sounds that triggered detection
DataSnipStruc = dir(fullfile(AudioDataPath, [DataFile(1:16) '*snippets/*.wav']));

if TranscExtract && nargin<2
    Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';
end
if TranscExtract && nargin<3
    % Set the path to logger data
    Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'logger',['20' Date]);
end

% Set the path to a working directory on the computer so logger data are
% transfered there and directly accessible for calculations
if TranscExtract
    WorkDir = ['~' filesep 'WorkingDirectory'];
end

%% Get the sample stamp of the detected vocalizations
fprintf(1,'*** Getting events for that day ***\n');
DataFileStruc = dir(fullfile(AudioDataPath, [DataFile(1:16) '*events.txt']));
Fid_Data = fopen(fullfile(DataFileStruc.folder,DataFileStruc.name));
EventsHeader = textscan(Fid_Data, '%s\t%s\t%s\t%s\t%s\t%s\t%s\n',1);
for hh=1:length(EventsHeader)
    if strfind(EventsHeader{hh}{1}, 'SampleStamp')
        EventsStampCol = hh;
    elseif strfind(EventsHeader{hh}{1}, 'Type')
        EventsEventTypeCol = hh;
    elseif strfind(EventsHeader{hh}{1}, 'FoodPortFront')
        EventsFoodPortFrontCol = hh;
    elseif strfind(EventsHeader{hh}{1}, 'FoodPortBack')
        EventsFoodPortBackCol = hh;
    elseif strfind(EventsHeader{hh}{1}, 'TimeStamp(s)')
        EventsTimeCol = hh;
    elseif strfind(EventsHeader{hh}{1}, 'Delay2Reward')
        EventsRewardCol = hh;
    end
end
Events = textscan(Fid_Data, '%s\t%f\t%s\t%s\t%f\t%f\t%f');
fclose(Fid_Data);
VocId = find(strcmp('Vocalization', Events{EventsEventTypeCol}));

%% Plot the cumulative number of triggers along time
fprintf(1,'*** Plotting cumulative events for that day ***\n');
ColorCode = get(groot,'DefaultAxesColorOrder');
ReTriggerVocId = intersect(VocId, find(~isnan(Events{EventsRewardCol})));
ReVocId = intersect(VocId, find(~(isnan(Events{EventsRewardCol}) + isinf(Events{EventsRewardCol}))));

F=figure(100);
plot(Events{EventsTimeCol}(VocId)/60,1:length(VocId), 'k-', 'Linewidth',2)
hold on
plot(Events{EventsTimeCol}(ReTriggerVocId)/60, 1:length(ReTriggerVocId), '-','Color',ColorCode(1,:),'Linewidth',2)
hold on
plot(Events{EventsTimeCol}(ReVocId)/60, 1:length(ReVocId), '-','Color',ColorCode(3,:),'Linewidth',2)
legend('Sound detection events', 'Sound trigger events', 'Rewarded sound events','Location','NorthWest')
xlabel('Time (min)')
ylabel('Cumulative sum of events')
hold off
title(sprintf('Subjects: %s  Date: %s  Time: %s', DataFile(1:4), DataFile(6:11), DataFile(13:16)))
hold off
saveas(F,fullfile(AudioDataPath,sprintf('%s_CumTrigger.fig', DataFile(1:16))))
saveas(F,fullfile(AudioDataPath,sprintf('%s_CumTrigger.jpeg', DataFile(1:16))))

%% This section is getting ready a plot around the time of one given detected vocalization
% % Calculate the offset of each soundfile output from the begining of the
% % task in number of samples. Offset_Y is the position of the first sample
% % of each file in the continuous recording (because we are dropping some
% % samples at each file change)
% Length_Filename = fullfile(AudioDataPath, sprintf('%s_Offset_Y.m',DataFile(1:16)));
% if ~exist(Length_Filename, 'file')
%     fprintf(1,'Calculating offset of each sound file to allign extract...\n')
%     FileChangeId = find(strcmp('ChangeFile', Events{EventsEventTypeCol}));
%     Nfiles = length(WavFileStruc);
%     Offset_Y = nan(Nfiles,1);
%     Length_Y = nan(Nfiles,1);
%     ExpDelay_Y = nan(Nfiles,1);
%     for yy=1:Nfiles
%         fprintf(1,'File %d/%d\n', yy,Nfiles)
%         % get the files in the correct order
%         Wavefile=dir(fullfile(WavFileStruc(yy).folder, sprintf('%s*_%d.wav',WavFileStruc(yy).name(1:(end-7)),yy)));
%         Wavefile_local = fullfile(WavFileStruc(yy).folder, Wavefile.name);
%         [Y,FS] = audioread(Wavefile_local);
%         if yy==1 % only calculate the length of the file, there is no offset
%             Length_Y(yy) = length(Y);
%             Offset_Y(yy) = 1;
%             ExpDelay_Y(yy) = 0;
%         else
%             % store the file length
%             Length_Y(yy) = length(Y);
%             % find the first snip of that sequence
%             for ss=1:length(DataSnipStruc)
%                 IndStamp1 = strfind(DataSnipStruc(ss).name, '_');
%                 IndStamp_last = IndStamp1(end);
%                 Sequence = str2double(DataSnipStruc(ss).name((IndStamp1(end-1)+1): (IndStamp1(end)-1)));
%                 Stamp=[];
%                 Ysnip = [];
%                 if Sequence == yy
%                     [Ysnip,~] = audioread(fullfile(DataSnipStruc(ss).folder, DataSnipStruc(ss).name));
%                     Stamp = str2double(DataSnipStruc(ss).name((IndStamp_last+1):end-3)); % Stamp is the position of the first sample of Ysnip
%                     break
%                 end
%             end
%             if Stamp<0
%                 Stamp = 2*2147483647 + Stamp;
%             end
%
%             if isempty(Stamp)
%             % Find the delay calculated by matlab when changing of file
%             ExpDelay_sec = Events{EventsRewardCol}(FileChangeId(yy-1));
%             ExpDelay_Y(yy) = round(ExpDelay_sec*FS);
%
%             % Align that snip with the wav of the sequence and get the stamp of the snip
%             Hyp_Stamp = Stamp-(sum(Length_Y(1:(yy-1))) + sum(Offset_Y(1:(yy-1))) + ExpDelay_Y(yy)); % This is the hypothetical position of the first sample of Ysnip in the present recording
%             Y_extract = Y((Hyp_Stamp-length(Ysnip)) : (Hyp_Stamp+2*length(Ysnip)));
%             Ncorr = length(Y_extract) - length(Ysnip);
%             Xcor = nan(Ncorr,1);
%             for cc=1:Ncorr
%                 Xcor(cc) = corr(Ysnip,Y_extract(cc:(cc+length(Ysnip)-1))); % Running a cross correlation between the raw signal centered around the hypothetize position of the snippet and the snippet of sounds
%             end
%             [~,I] = max(abs(Xcor)); % I is the first index of Y_extract that is best alligned with the first index of Ysnip
%             True_Stamp = Hyp_Stamp-length(Ysnip) + I; % This is the index of the first element of Ysnip in Y
%             Offset_Y(yy) = Stamp-True_Stamp+1; % This is the position of the first element of Y in the continuous recording
%             figure(2)
%             cla
%             subplot(1,2,1)
%             plot(1:Ncorr,Xcor)
%             hold on
%             line([I I], [min(Xcor) max(Xcor)], 'Color','r', 'LineStyle','--');
%             text(I , mean(Xcor), sprintf('%d',I))
%             xlabel('Lag')
%             ylabel('correlation coefficient')
%             hold off
%             subplot(1,2,2)
%             plot(Y((Stamp-Offset_Y(yy)) : (Stamp-Offset_Y(yy)+ length(Ysnip)-1)), 'k-', 'LineWidth',2)
%             hold on
%             plot(Ysnip, 'r-')
%             hold off
%             legend('Continuous recording','SoundSnippet')
%             title(sprintf('File %d/%d\n', yy,Nfiles))
%             pause()
%         end
%     end
%     save(Length_Filename,'Offset_Y', 'Length_Y')
% end

% Plot the results
% List of all the detected calls
% Chose a vocalization to center the vizualization tool
if PlotIndivFile
    % first get the length of all files recorded to localize extract within the
    % files
    Subj= DataFile(1:4);
    Time = DataFile(13:16);
    [Length_Y] = get_raw_file_length(AudioDataPath, Subj, Date, Time);
    
    fprintf(1,'Now plot the results around a given call\n')
    if ~exist('FS', 'var')
        [~,FS] = audioread(fullfile(WavFileStruc(1).folder, WavFileStruc(1).name));
    end
    Buffer = 1*FS;% We choose to have one second of recording before and after sound detection onset
    IndCenterVoc=1;
    while ~isempty(IndCenterVoc)
        FullStamps = Events{EventsStampCol}(VocId);
        
        SoundtimeMinSec = cell(length(FullStamps),1);
        fprintf(1,'The following vocalizations where automatically detected by voc operant.\n Chose the index of the vocalization you want to look at:\n')
        for ss=1:length(FullStamps)
            Ind_ = strfind(FullStamps{ss}, '_');
            Stamp = str2double(FullStamps{ss}((Ind_+1):end));
            if Stamp<0
                Stamp = 2*2147483647 + Stamp; % Correcion of soundmexpro bug that coded numbers in 32 bits instead of 64bits
            end
            MinStamp = floor(Stamp/(60*FS));
            SecStamp = (Stamp - (MinStamp*60*FS))/FS;
            SoundtimeMinSec{ss} = sprintf('%dmin %.1fs', MinStamp, SecStamp);
            fprintf(1, '%d. %s\n', ss,SoundtimeMinSec{ss})
        end
        IndCenterVoc = input('Index of your choice (leave empty to quit):\n');
        if isempty(IndCenterVoc)
            break
        end
        
        Ind_ = strfind(FullStamps{IndCenterVoc}, '_');
        Seq = str2double(FullStamps{IndCenterVoc}(1:(Ind_-1)));
        Stamp = str2double(FullStamps{IndCenterVoc}((Ind_+1):end));
        if Stamp<0
            Stamp = 2*2147483647 + Stamp; % Correcion of soundmexpro bug that coded numbers in 32 bits instead of 64bits
        end
        
        % Get the recording data for the corresponding sequence
        WavFileStruc_local = dir(fullfile(AudioDataPath, sprintf('%s*mic*_%d.wav',DataFile(1:16), Seq)));
        try
            Wavefile_local = fullfile(WavFileStruc_local.folder, WavFileStruc_local.name);
            [Y,FS] = audioread(Wavefile_local);
        catch
            fprintf(1,'Warning: the audiofile %s cannot be read properly and will not be plotted\n', Wavefile_local);
            Y = 0;
        end
        Y_section_beg = max(1,Stamp - sum(Length_Y(1:(Seq-1))) - Buffer); % Make sure we don't request before the beginning of the raw wave file
        Pre_stamp = min(Buffer, Stamp - sum(Length_Y(1:(Seq-1)))); % Length of the sound section before sound detection onset
        Y_section_end = min(length(Y), Stamp - sum(Length_Y(1:(Seq-1))) + Buffer); % Make sure we don't request after the end of the aw wave file
        Post_stamp = min(Buffer, length(Y)- (Stamp - sum(Length_Y(1:(Seq-1)))));% Length of the sound section after sound detection onset
        Y_section = Y(Y_section_beg:Y_section_end);
        
        % Plot the waveforms of the recording around the stamp of the vocalization
        figure(2)
        cla
        plot(Y_section, 'Color', 'k')
        
        
        % plot the sound snippets in red on top
        % Get the sound snippets from the sounds that triggered detection
        DataSnipStruc = dir(fullfile(AudioDataPath, sprintf('%s*snippets/*snipfile_%d_*.wav', DataFile(1:16), Seq)));
        NbSnip_local = length(DataSnipStruc);
        SnipDur = nan(NbSnip_local,1);
        Stamp_all = nan(NbSnip_local,1);
        VocId_all = nan(NbSnip_local,1);
        Voc_event_nb = 0;
        for ss=1:NbSnip_local
            IndStamp1 = strfind(DataSnipStruc(ss).name, '_');
            IndStamp_last = IndStamp1(end);
            Stamp_local = str2double(DataSnipStruc(ss).name((IndStamp_last+1):end-4));
            
            if (Stamp_local < (Stamp + Post_stamp)) && (Stamp_local > (Stamp - Pre_stamp))
                [Ysnip,FS] = audioread(fullfile(DataSnipStruc(ss).folder, DataSnipStruc(ss).name));
                %             Stamp_local_context = Stamp_local - (Stamp - Pre_stamp)+1;% This is where the first sample of the sound Ysnip should be in Y_section
                % There are often lay-off between the sample value and the
                % actual position within the recording, estimating that lay-off
                % using cross correlation
                DiffY = length(Y_section)-length(Ysnip);
                XcorrY=nan(1,DiffY+1);
                for cc=0:DiffY
                    XcorrY(cc+1) = transpose(Y_section(cc+(1:length(Ysnip)))) * Ysnip;
                end
                [~,Stamp_local_context] = max(abs(XcorrY));
                Sequence = str2double(DataSnipStruc(ss).name((IndStamp1(end-1)+1):(IndStamp1(end)-1)));
                Voc_event_nb = Voc_event_nb + 1;
                figure(2)
                hold on
                plot(Stamp_local_context:(Stamp_local_context+length(Ysnip)-1), Ysnip, 'Color', 'g--')
                % Search for the corresponding datapoint in the event log
                Index =[];
                for ii=1:length(VocId)
                    if strcmp(Events{EventsStampCol}{VocId(ii)}, sprintf('%d_%d', Sequence, Stamp_local))
                        Index = ii;
                        break
                    end
                end
                if isempty(Index)
                    error('event %d_%d cannot be found in the eventlog file\n', Sequence, Stamp_local);
                else
                    Stamp_all(Voc_event_nb) = Stamp_local;
                    SnipDur(Voc_event_nb) = length(Ysnip);
                    VocId_all(Voc_event_nb) = VocId(Index);
                end
            end
        end
        
        % plot the time stamp when sound was detected (last time point of each sound snippet)
        % color code the point given the reward status
        % Get the VocId should be plotted, from the log
        VocId_log = nan(1,length(VocId));
        Stamp_all_log = nan(1,length(VocId));
        SnipDur_log = nan(1,length(VocId));
        Voc_event_nb_log = 0;
        for ii=1:length(VocId)
            FullStamp = Events{EventsStampCol}{VocId(ii)};
            IndStamp = strfind(FullStamp, '_');
            Seq_local = str2double(FullStamp(1:(IndStamp-1)));
            Stamp_local = str2double(FullStamp((IndStamp+1):end));
            if (Seq_local == Seq) && (Stamp_local >= (Stamp-Pre_stamp)) && (Stamp_local <= (Stamp+Post_stamp)) % This is a vocalization from the same section
                Voc_event_nb_log = Voc_event_nb_log+1;
                VocId_log(Voc_event_nb_log) = VocId(ii);
                Stamp_all_log(Voc_event_nb_log) = Stamp_local;
                try
                    SnipDur_log(Voc_event_nb_log) = SnipDur(find(VocId_all==VocId(ii)));
                catch
                    warning('This sound detection event does not have a corresponding sound snippet %s',FullStamp);
                    SnipDur_log(Voc_event_nb_log) = NaN;
                end
            end
        end
        VocId_log = VocId_log(1:Voc_event_nb_log);
        Stamp_all_log = Stamp_all_log(1:Voc_event_nb_log);
        SnipDur_log = SnipDur_log(1:Voc_event_nb_log);
        
        ReFront = Events{EventsFoodPortFrontCol}(VocId_log);
        ReBack = Events{EventsFoodPortBackCol}(VocId_log);
        Colorcode = get(groot,'DefaultAxesColorOrder');
        ColorFront= nan(size(ReFront,1),3);
        ColorBack= ColorFront;
        ColorFront(find(ReFront==1),:) = repmat(Colorcode(3,:),sum(ReFront==1),1);
        ColorBack(find(ReBack==1),:) = repmat(Colorcode(3,:), sum(ReBack==1),1);
        ColorFront(find(ReFront==0),:) = repmat([0 0 0],sum(ReFront==0),1);
        ColorBack(find(ReBack==0),:) = repmat([0 0 0], sum(ReBack==0),1);
        ColorFront(find(isnan(ReFront)),:) = repmat([1 0 0],sum(isnan(ReFront)),1);
        ColorBack(find(isnan(ReBack)),:) = repmat([1 0 0], sum(isnan(ReBack)),1);
        figure(2)
        hold on
        sc=scatter(Stamp_all_log+SnipDur_log-1 - (Stamp - Pre_stamp), 1.2*ones(length(VocId_log),1), 10, ColorFront, 'filled');
        hold on
        sc=scatter(Stamp_all_log+SnipDur_log-1 - (Stamp - Pre_stamp), 1.6*ones(length(VocId_log),1), 10, ColorBack, 'filled');
        
        XaxisRes=20;
        LengthY = length(Y_section);
        set(gca, 'XLim',[0;LengthY]);
        
        if (LengthY/FS)>XaxisRes
            XaxisStep = round(LengthY/(FS*XaxisRes));
        else
            XaxisStep = round(LengthY/(FS*XaxisRes),2);
        end
        set(gca,'XTick', 0:(XaxisStep*FS):LengthY, 'XTickLabel', 0:XaxisStep:round(LengthY/FS), 'YLim',[-1 2]);
        xlabel('Time (s)');
        ylabel('Sound pressure')
    end
    
end

%% Extracting sound events
% The samplestamp given by sound mex is not really reliable, so for each
% sound snippet, you want to find its exact location in the continuous
% recording files, then using TTL pulses, retrieve the time it correspond
% to in Deuteron, if requested.

% Checking what we have in terms of vocalization localization/extraction
ExpStartTime = DataFile(13:16);
VocExt_dir = dir(fullfile(AudioDataPath,sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)));

% Then run the logger extraction, allignment, and vocalization extraction
if TranscExtract
    fprintf(1,'*** Extract Logger data if not already done ***\n');
    % Find the ID of the recorded bats
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:P200','basic');
    RowData = find((cell2mat(RecTableData(2:end,1))== str2double(Date))) +1;
    DataInfo = RecTableData(RowData,:);
    Header = RecTableData(1,:);
    BatIDCol = find(contains(Header, 'Bat'));
    
    % extract logger data if not already done
    All_loggers_dir = dir(fullfile(Logger_dir, '*ogger*'));
    DirFlags = [All_loggers_dir.isdir];
    % Extract only those that are directories.
    All_loggers_dir = All_loggers_dir(DirFlags);
    TransceiverReset = struct(); % These are possible parameters for dealing with change of transceiver or sudden transceiver clock change. Set to empty before the first extraction
    LoggerName = cell(length(All_loggers_dir),1);
    BatID = cell(length(All_loggers_dir),1);
    for ll=1:length(All_loggers_dir)
        Logger_i = fullfile(Logger_dir,All_loggers_dir(ll).name);
        Ind = strfind(All_loggers_dir(ll).name, 'r');
        Logger_num = str2double(All_loggers_dir(ll).name((Ind+1):end));
        NLogCol = find(contains(Header, 'NL'));% Columns of the neural loggers
        ALogCol = find(contains(Header, 'AL'));% Columns of the audio loggers
        LogCol = NLogCol(find(cell2mat(DataInfo(NLogCol))==Logger_num));
        if isempty(LogCol) % This is an audiologger and not a neural logger
            LogCol = ALogCol(find(cell2mat(DataInfo(ALogCol))==Logger_num));
            LoggerName{ll} = ['AL' num2str(Logger_num)];
        else
            LoggerName{ll} = ['NL' num2str(Logger_num)];
        end
        BatID{ll} = DataInfo{BatIDCol(find(BatIDCol<LogCol,1,'last'))};
        ParamFiles = dir(fullfile(Logger_i,'extracted_data','*extract_logger_data_parameters*mat'));
        if isempty(ParamFiles) || ForceExtract
            fprintf(1,'-> Extracting %s\n',All_loggers_dir(ll).name);
            
            % Bring data back on the computer
            Logger_local = fullfile(WorkDir, All_loggers_dir(ll).name);
            fprintf(1,'Transferring data from the server %s\n on the local computer %s\n', Logger_i, Logger_local);
            mkdir(Logger_local)
            [s,m,e]=copyfile(Logger_i, Logger_local, 'f');
            if ~s
                m %#ok<NOPRT>
                e %#ok<NOPRT>
                error('File transfer did not occur correctly for %s\n', Logger_i);
            end
            
            % run extraction
            if Logger_num==16 && str2double(Date)<190501
                % extract_logger_data(Logger_local, 'BatID', num2str(BatID), 'ActiveChannels', [0 1 2 3 4 5 6 7 8 9 10 12 13 14 15], 'AutoSpikeThreshFactor',5,'TransceiverReset',TransceiverReset)
                extract_logger_data(Logger_local, 'BatID', num2str(BatID{ll}), 'ActiveChannels', [0 1 2 3 4 5 6 7 8 9 10 12 13 14 15],'TransceiverReset',TransceiverReset,'AutoSpikeThreshFactor',4)
            else
                %extract_logger_data(Logger_local, 'BatID', num2str(BatID),'TransceiverReset',TransceiverReset)
                extract_logger_data(Logger_local, 'BatID', num2str(BatID{ll}),'TransceiverReset',TransceiverReset,'AutoSpikeThreshFactor',4)
            end
            
            % Keeps value of eventual clock reset
            Filename=fullfile(Logger_local, 'extracted_data', sprintf('%s_20%s_EVENTS.mat', num2str(BatID{ll}),Date));
            NewTR = load(Filename, 'TransceiverReset');
            if ~isempty(fieldnames(NewTR.TransceiverReset))% this will be used in the next loop!
                TransceiverReset = NewTR.TransceiverReset;
            end
            
            % Bring back data on the server
            fprintf(1,'Transferring data from the local computer %s\n back on the server %s\n', Logger_i, Logger_local);
            Remote_dir = fullfile(Logger_i, 'extracted_data');
            mkdir(Remote_dir)
            [s,m,e]=copyfile(fullfile(Logger_local, 'extracted_data'), Remote_dir, 'f');
            if ~s
                TicTransfer = tic;
                while toc(TicTransfer)<30*60
                    [s,m,e]=copyfile(fullfile(Logger_local, 'extracted_data'), Remote_dir, 'f');
                    if s
                        return
                    end
                end
                if ~s
                    s %#ok<NOPRT>
                    m %#ok<NOPRT>
                    e %#ok<NOPRT>
                    error('File transfer did not occur correctly for %s\n Although we tried for 30min\n', Remote_dir);
                else
                    fprintf('Extracted data transfered back on server in:\n%s\n',  Remote_dir);
                end
            else
                fprintf('Extracted data transfered back on server in:\n%s\n',  Remote_dir);
            end
            if s  %erase local data
                [sdel,mdel,edel]=rmdir(WorkDir, 's');
                if ~sdel
                    TicErase = tic;
                    while toc(TicErase)<30*60
                        [sdel,mdel,edel]=rmdir(WorkDir, 's');
                        if sdel
                            return
                        end
                    end
                end
                if ~sdel
                    sdel %#ok<NOPRT>
                    mdel %#ok<NOPRT>
                    edel %#ok<NOPRT>
                    error('File erase did not occur correctly for %s\n Although we tried for 30min\n', WorkDir);
                end
            end
            
        else
            fprintf(1,'-> Already done for %s\n',All_loggers_dir(ll).name);
        end
    end
    
    % Get the serial numbers of the audiologgers that the two implanted
    % bats wear
    NLCol = find(contains(Header, 'NL'));
    ALThroatCol = find(contains(Header, 'AL-throat'));
    SerialNumberAL = nan(length(NLCol),1);
    SerialNumberNL = nan(length(NLCol),1);
    for dd=1:length(NLCol)
        SerialNumberAL(dd) = DataInfo{ALThroatCol(find(ALThroatCol<NLCol(dd),1,'last'))};
        SerialNumberNL(dd) = DataInfo{NLCol(dd)};
    end
    
    % Alligning TTL pulses between soundmexpro and Deuteron
    % for the Operant session
    
    TTL_dir = dir(fullfile(AudioDataPath,sprintf( '%s_%s_TTLPulseTimes.mat', Date, ExpStartTime)));
    if isempty(TTL_dir) || ForceAllign
        fprintf(1,'*** Alligning TTL pulses for the operant session ***\n');
        align_soundmexAudio_2_logger(AudioDataPath, Logger_dir, ExpStartTime,'TTL_pulse_generator','Avisoft','Method','risefall', 'Session_strings', {'all voc reward start', 'all voc reward stop'}, 'Logger_list', [SerialNumberAL; SerialNumberNL]);
    else
        fprintf(1,'*** ALREADY DONE: Alligning TTL pulses for the operant session ***\n');
    end
    if isempty(VocExt_dir) || ForceVocExt1
        fprintf(1,'*** Localizing and extracting vocalizations that triggered the sound detection ***\n');
        voc_localize_operant(AudioDataPath, DataFile(1:4),Date, ExpStartTime, 'UseSnip',0)
    else
        fprintf(1,'*** ALREADY DONE: Localizing and extracting vocalizations that triggered the sound detection ***\n');
    end
    
    load(fullfile(AudioDataPath, sprintf('%s_%s_VocExtractTimes.mat', Date, ExpStartTime)),'Re_transc_time');
    if ~exist('Re_transc_time','var') || ForceVocExt3
        fprintf(1,'*** Localizing and extracting reward times ***\n');
        voc_localize_operant_patch(AudioDataPath, DataFile(1:4),Date, ExpStartTime)
    end
    %% Identify the same vocalizations on the piezos and save sound extracts, onset and offset times
    fprintf(' LOCALIZING VOCALIZATIONS ON PIEZO RECORDINGS\n')
    LogVoc_dir = dir(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)));
    if isempty(LogVoc_dir) || ForceVocExt1 || ForceVocExt2 || ForceVocExt3
        get_logger_data_voc(AudioDataPath, Logger_dir,Date, ExpStartTime, 'SerialNumber',SerialNumberAL);
    else
        fprintf(1,'Using already processed data\n')
        
    end
    
    %% Identify who is calling
    fprintf(' IDENTIFY WHO IS CALLING\n')
    WhoCall_dir = dir(fullfile(Logger_dir, sprintf('*%s_%s*whocalls*', Date, ExpStartTime)));
    if isempty(WhoCall_dir) || ForceVocExt1 || ForceWhoID || ForceVocExt2
        who_calls_patch(AudioDataPath, Logger_dir,Date, ExpStartTime,200,1,1,0,'Working_dir','/Users/elie/Documents/WorkingWhoPatch');
    else
        fprintf(1,'Using already processed data\n')
    end
    % Save the ID of the bat for each logger
    save(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, 200)), 'BatID','LoggerName','-append')

     %% Explore what is said
    fprintf('\n*** Identify what is said ***\n')
    WhatCall_dir = dir(fullfile(Logger_dir,'VocExtracts', sprintf('*%s_%s*Elmt*Raw.wav', Date, ExpStartTime)));
    if isempty(WhatCall_dir) || ForceVocExt1 || ForceVocExt2 || ForceWhat
        what_calls(Logger_dir,Date, ExpStartTime);
    else
        fprintf('\n*** ALREADY DONE: Identify what is said ***\n')
    end
    
elseif isempty(VocExt_dir) || ForceVocExt1
    fprintf(1,'*** Localizing and extracting vocalizations that triggered the sound detection ***\n');
    fprintf(1,'NOTE: no transceiver time extraction\n')
    voc_localize_operant(AudioDataPath, DataFile(1:4),Date, ExpStartTime, 'UseSnip',0,'TransceiverTime',0)
else
    fprintf(1,'*** ALREADY DONE: Localizing and extracting vocalizations that triggered the sound detection ***\n');
end





end
