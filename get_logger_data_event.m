function [AllActions, UActionText]= get_logger_data_event(ExtData_dir, Date, Buffer, LoggerID, Output_dir)
% AllActions is a cell array of length the number of actions containing
% each 2 columns giving the onset and offset times in ms of behavioral
% events
if nargin<3
    Buffer = 500; %ms
end

AllSep = strfind(ExtData_dir, filesep);
if nargin<4
    if length(ExtData_dir)==AllSep(end)
        LoggerID = ExtData_dir((AllSep(end-2)+1):(AllSep(end-1)-1));
    else
        LoggerID = ExtData_dir((AllSep(end-1)+1):(AllSep(end)-1));
    end
end
if nargin<5
    if length(ExtData_dir)==AllSep(end)
        Output_dir = ExtData_dir(1:(AllSep(end-2)));
    else
        Output_dir = ExtData_dir(1:(AllSep(end-1)));
    end
end
fprintf(1, 'Figures and data are saved to %s\n', Output_dir)
MaxEventDur = 1000; %ms
EventfilePath = dir(fullfile(ExtData_dir, sprintf('*%s*EVENTS.mat', Date)));
Eventfile = load(fullfile(EventfilePath.folder, EventfilePath.name));

%% identify all the free text strings reported
FreeTextInd = find(contains(Eventfile.event_types_and_details, 'Free text'));
UFreeText_all = unique(Eventfile.event_types_and_details(FreeTextInd)); %#ok<FNDSB> % unique labels
% Remove 'Free Text' and isolate actions items
SpaceInd = strfind(UFreeText_all, ' ');
StartInd= strfind(UFreeText_all, 'start');
StopInd= strfind(UFreeText_all, 'stop');
UActionText = cell(size(UFreeText_all));
UFullActionText = cell(size(UFreeText_all));
for ii=1:length(SpaceInd)
    UFullActionText{ii} = UFreeText_all{ii}((SpaceInd{ii}(2)+1) :end);
    if ~isempty(StartInd{ii})
        UActionText{ii} = UFreeText_all{ii}((SpaceInd{ii}(2)+1) :(StartInd{ii}-1));
    elseif ~isempty(StopInd{ii})
        UActionText{ii} = UFreeText_all{ii}((SpaceInd{ii}(2)+1) :(StopInd{ii}-1));
    else
        UActionText{ii} = UFreeText_all{ii}((SpaceInd{ii}(2)+1):end);
    end
end
UActionText = unique(UActionText);

%% construct an ethogram of the recording
NAction = length(UActionText);
AllActions = cell(NAction,1);
% Find the reference time to plot the ethogram
RefTime = Eventfile.event_timestamps_usec(find(contains(Eventfile.event_types_and_details, 'Started recording'),1))*10^-3;
% Get the indices of all the starting points
IndStart = find(contains(Eventfile.event_types_and_details, 'start'));
% Get the indices of all the stoping points
IndStop = find(contains(Eventfile.event_types_and_details, 'stop'));
Fig=figure();
ColorCode = get(groot,'DefaultAxesColorOrder');
ColorCode = [ColorCode(1:6,:) ; 0 0 1; ColorCode(7,:); 0.85 0.6940 0.556; 0 0 0; 1 0 0; 0.301 0.078 0.741; ];
% loop through each action
for aa=1:NAction
    IndAction = find(contains(Eventfile.event_types_and_details, UActionText{aa}));
    IndActionStart = intersect(IndAction, IndStart);
    if isempty(IndActionStart) % This action is ponctual
        AllActions{aa} = Eventfile.event_timestamps_usec(IndAction)*10^-3;
        hold on
        plot((AllActions{aa}-RefTime)*10^-3, aa, '*','Color', ColorCode(aa,:), 'MarkerSize', 5)
    else
        IndActionStop = intersect(IndAction, IndStop);
        AllActions{aa} = nan(max(length(IndActionStop), length(IndActionStart)),2);
        % loop through action starts and figure out if a stop is closer to
        % it than the next start
        e_count = 0;
        for ii=1:length(IndActionStart)
            % get the stop closer to the start
            mm = find(min(IndActionStop - IndActionStart(ii)));
            if (IndActionStop(mm) - IndActionStart(ii)) <0
                error('Error in parsing the action, the stop happens before the start\n')
            end
            if (ii==length(IndActionStart)) || (IndActionStop(mm) - IndActionStart(ii)) < diff(IndActionStart(ii:(ii+1)))  % this stop happens before another start, all is good save that pair of event
                e_count= e_count+1;
                AllActions{aa}(e_count,1) = Eventfile.event_timestamps_usec(IndActionStart(ii))*10^-3;
                AllActions{aa}(e_count,2) = Eventfile.event_timestamps_usec(IndActionStop(mm))*10^-3;
                % plot the action
                hold on
                plot((AllActions{aa}(e_count,:)-RefTime)*10^-3, aa*ones(2,1), '-', 'Color', ColorCode(aa,:), 'LineWidth',2);
                if ii<length(IndActionStart)
                    % Find the next start
                    nn = find(min(IndActionStart(ii+1:end) - IndActionStop(mm)));
                    if (IndActionStart(ii+nn) - IndActionStop(mm)) < diff(IndActionStop(mm:(mm+1)))
                        % The next start happens before another stop, all is
                        % good, proceed to the next start and supress that stop
                        IndActionStop(mm) = [];
                    else % There are several stops in a row, get rid of them up to the next start
                        SupInd = find(diff(IndActionStop(mm:end))<(IndActionStart(ii+nn) - IndActionStop(mm)));
                        IndActionStop(mm : (mm+ SupInd)) = [];
                    end
                end
            else 
                % This start misses a stop, do not consider it consider the next start
            end
        end
        AllActions{aa} = AllActions{aa}(1:e_count,:);
    end
end

set(gca,'YLim',[0 NAction+1], 'YTick',1:NAction, 'YTickLabel', UActionText)
xlabel('Time (in s)')
title(sprintf('Ethogram Implanted bat %s', Date))

orient(Fig,'landscape')
Fig.PaperPositionMode = 'auto';
set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
set(Fig,'PaperOrientation','landscape');
set(Fig,'PaperUnits','normalized');
set(Fig,'PaperPosition', [0 0 1 1]);

print(Fig,fullfile(Output_dir,sprintf('%s_%s_Actuogram.pdf', Date, LoggerID)),'-dpdf')

%% extract the neural activity during the whole recording time
PlayBackOffset = AllActions{contains(UActionText, 'playback')}(end,2);
fprintf('Extract Single Unit spike arrival times for the whole experiment session\n')
Delay = 0;% time in ms to extract data before and after event onset/offset
[~,~, ~, SpikeSU] = extract_timeslot_LFP_spikes(ExtData_dir, [RefTime, PlayBackOffset], Delay, NaN,[0 0 0 1]);
save(fullfile(Output_dir, sprintf('%s_SpikeSU_RecordPlayback.mat', Date)),'SpikeSU');
Response_samprate = 1;% Sampling rate of the KDE in Hz
Bin_ms = 1000; % size of the KDE binning
t=-Delay: Bin_ms : round((PlayBackOffset-RefTime + Delay)/Bin_ms)*Bin_ms;
Psth_KDEfiltered_SpikeSU = cell(1,length(SpikeSU));
Psth_KDEfiltered_SpikeSU_t = cell(1,length(SpikeSU));
for uu=1:length(SpikeSU)
    fprintf('calculate spike density estimate SU %d/%d\n', uu, length(SpikeSU))
    % calculate the density estimate (Kernel density estimation: L2 risk
    % minimization, with reflection boundary)
    [y,Psth_KDEfiltered_SpikeSU_t{uu},~]=ssvkernel(SpikeSU{uu},t);
%     [Psth_KDEfiltered_SpikeSU_t{uu},y]=Bayes_v1(SpikeSU{uu});
    % y is a density function that sums to 1
    % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 1 s)
    % multiplying by the response sampling rate in Hz gives the expected spike rate in spike/s
    Psth_KDEfiltered_SpikeSU{uu} =  y * length(SpikeSU{uu}) * Response_samprate;
end

Fig=figure();
% Plot KDE
SpikeCol = [1 0.2 0.2 0.3;  0.5 0.7 0.5 0.8; 0.9 0.1 0.2 0.3; [get(groot,'DefaultAxesColorOrder') repmat(0.3,7,1)]];
LegendKDE ={};
for uu=1:length(SpikeSU)
    yyaxis right
    hold on
    plot(Psth_KDEfiltered_SpikeSU_t{uu}*10^-3,zscore(Psth_KDEfiltered_SpikeSU{uu}),'-','Color',SpikeCol(uu,:), 'LineWidth',2, 'DisplayName', sprintf('Single Unit %d',uu))
    LegendKDE{uu}  = sprintf('Single Unit %d',uu);
    hold off    
end
ylabel('Zscored Spike Rate (/s)')
legend(LegendKDE, 'AutoUpdate','off')
YlimRight = get(gca, 'YLim');
set(gca, 'YLim', [0 YlimRight(2)*2])

% Plot actions on top
yyaxis left
for aa=1:NAction
    if size(AllActions{aa},2)==1
        hold on
        plot((AllActions{aa}-RefTime)*10^-3, aa, '*','Color', ColorCode(aa,:), 'MarkerSize', 5)
    elseif size(AllActions{aa},2)==2
        for e_count=1:size(AllActions{aa},1)
            hold on
            plot((AllActions{aa}(e_count,:)-RefTime)*10^-3, aa*ones(2,1), '-', 'Color', ColorCode(aa,:), 'LineWidth',2);
        end
    end
end
xlabel('Time (in s)')
title(sprintf('Ethogram and SU bat %s', Date))

% Plot spikes!
LegendPlot = cell(length(SpikeSU) + length(UActionText),1);
LegendPlot((length(SpikeSU)+1) : end) = UActionText;
yyaxis left
for uu=1:length(SpikeSU)
    hold on
    plot(SpikeSU{uu}*10^-3, ones(length(SpikeSU{uu}),1)*(-uu/(length(SpikeSU)+1)), '.', 'MarkerSize',5, 'Color', SpikeCol(uu,(1:3)))
    LegendPlot{length(SpikeSU)+1-uu} = sprintf('Single Unit %d',uu);
end
set(gca,'YLim',[-1 1]*(NAction+1), 'YTick',[-flip(1:length(SpikeSU))/(length(SpikeSU)+1) 1:NAction], 'YTickLabel', LegendPlot)

orient(Fig,'landscape')
Fig.PaperPositionMode = 'auto';
set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
set(Fig,'PaperOrientation','landscape');
set(Fig,'PaperUnits','normalized');
set(Fig,'PaperPosition', [0 0 1 1]);

print(Fig,fullfile(Output_dir,sprintf('%s_%s_NeuroActuogram.pdf', Date, LoggerID)),'-dpdf')

%% Extract the neural activity around each behavioral event
% Only work on a subset of the events
Response_samprate = 100;% Sampling rate of the KDE in Hz
Bin_ms = 10; % size of the KDE binning
Cleared_events = {'eating', 'female', 'playback','record only', 'food in', 'vocalization', 'yawning'};
Cleared_events_Ind = nan(length(Cleared_events),1);
for ee=1:length(Cleared_events)
    Cleared_events_Ind(ee) = find(contains(UActionText, Cleared_events{ee}));
end
GoodEvents = setdiff(1:length(UActionText), Cleared_events_Ind);
NAction = length(GoodEvents);
SpikesTimesSU_Action = cell(NAction,1);
SpikesTimesT_Action = cell(NAction,1);
LFP_Action =  cell(NAction,1);
% Raw_Action = cell(NAction,1);
Event_duration =  cell(NAction,1);
Psth_KDEfilteredT_Action_scalef = cell(NAction,1);
Psth_KDEfilteredT_Action_t = cell(NAction,1);
Psth_KDEfilteredT_Action = cell(NAction,1);
Psth_KDEfilteredSU_Action_scalef = cell(NAction,1);
Psth_KDEfilteredSU_Action_t = cell(NAction,1);
Psth_KDEfilteredSU_Action = cell(NAction,1);
Average_Psth_KDEfiltered_TAction = cell(NAction,1);
Average_Psth_KDEfiltered_SUAction = cell(NAction,1);
Overall_Psth_KDEfilteredT_Action_t = cell(NAction,1);
Overall_Psth_KDEfilteredT_Action = cell(NAction,1);
 Overall_Psth_KDEfilteredT_Action_confb95 = cell(NAction,1);
Overall_Psth_KDEfilteredSU_Action_t = cell(NAction,1);
Overall_Psth_KDEfilteredSU_Action = cell(NAction,1);
 Overall_Psth_KDEfilteredSU_Action_confb95 = cell(NAction,1);
Overall_KDE_duration = nan(NAction,1);
for nn=1:NAction
    aa = GoodEvents(nn);
    % extracting neural data around the event
    fprintf(1,'Extracting neural data for %s events\n', UActionText{aa});
    [~ ,LFP_Action{nn}, SpikesTimesT_Action{nn},SpikesTimesSU_Action{nn}] = extract_timeslot_LFP_spikes(ExtData_dir, AllActions{aa}, Buffer, MaxEventDur, [0 1 1 1]);
    
    % calculate event duration
    Event_duration{nn} = AllActions{aa}(:,2) - AllActions{aa}(:,1);
    Event_duration{nn}(Event_duration{nn}>= (MaxEventDur + Buffer)) = MaxEventDur + Buffer;
    
    % calculate event PSTH
    % calculate the density estimate for each event and each unit
    Nevents = length(Event_duration{nn});
    NSU = size(SpikesTimesSU_Action{nn},2);
    NT = size(SpikesTimesT_Action{nn},2);
    Psth_KDEfilteredT_Action_t{nn} = cell(Nevents, NT);
    Psth_KDEfilteredT_Action{nn} = cell(Nevents, NT);
    Psth_KDEfilteredT_Action_scalef{nn} = nan(Nevents, NT);
    for ee=1:Nevents
        % Calculate the t for KDE
        t=-Buffer: Bin_ms : round((Event_duration{nn}(ee) + Buffer)/Bin_ms)*Bin_ms;
        % Calculate the kde for each Unit
        for uu=1:NT
            [y,Psth_KDEfilteredT_Action_t{nn}{ee,uu},~]=ssvkernel(SpikesTimesT_Action{nn}{ee,uu},t);
            % y is a density function that sums to 1
            % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 10 ms)
            % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
            Psth_KDEfilteredT_Action{nn}{ee,uu} =  y * length(SpikesTimesT_Action{nn}{ee,uu}) * Response_samprate/1000;
            Psth_KDEfilteredT_Action_scalef{nn}(ee,uu) = max(Psth_KDEfilteredT_Action{nn}{ee,uu});
        end
    end

    Psth_KDEfilteredSU_Action_t{nn} = cell(Nevents, NSU);
    Psth_KDEfilteredSU_Action{nn} = cell(Nevents, NSU);
    Psth_KDEfilteredSU_Action_scalef{nn} = nan(Nevents, NSU);
    for ee=1:Nevents
        % Calculate the t for KDE
        t=-Buffer: Bin_ms : round((Event_duration{nn}(ee) + Buffer)/Bin_ms)*Bin_ms;
        % Calculate the kde for each Unit
        for uu=1:NSU
            [y,Psth_KDEfilteredSU_Action_t{nn}{ee,uu},~]=ssvkernel(SpikesTimesSU_Action{nn}{ee,uu},t);
            % y is a density function that sums to 1
            % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 10 ms)
            % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
            Psth_KDEfilteredSU_Action{nn}{ee,uu} =  y * length(SpikesTimesSU_Action{nn}{ee,uu}) * Response_samprate/1000;
            Psth_KDEfilteredSU_Action_scalef{nn}(ee,uu) = max(Psth_KDEfilteredSU_Action{nn}{ee,uu});
        end
    end
    
    % calculate the density estimate across all events for each unit for
    % the common event duration
    Overall_KDE_duration(nn) = min(Event_duration{nn})+Buffer;
    t=-Buffer: Bin_ms : round(Overall_KDE_duration(nn)/Bin_ms)*Bin_ms;
    Overall_Psth_KDEfilteredT_Action_t{nn} = cell(NT,1);
    Overall_Psth_KDEfilteredT_Action_confb95{nn} = cell(NT,1);
    Overall_Psth_KDEfilteredT_Action{nn} = cell(NT,1);
    for uu=1:NT
            % the output of ssvkernel is a density function that sums to 1
            % multiplying by the total number of spikes gives the number of expected spike per time bin (here 10 ms)
            % dividing by the number of events gives the number of expected spike per time bin (here 10 ms) per event 
            % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
            try
                LocalSpikes = cell2mat(SpikesTimesT_Action{nn}(:,uu));
            catch
                LocalSpikes = cell2mat(SpikesTimesT_Action{nn}(:,uu)');
            end
            [y,Overall_Psth_KDEfilteredT_Action_t{nn}{uu},~, ~,~,Overall_Psth_KDEfilteredT_Action_confb95{nn}{uu}, ~]=ssvkernel(LocalSpikes,t);
            Overall_Psth_KDEfilteredT_Action{nn}{uu} =  y * length(LocalSpikes)/Nevents * Response_samprate/1000;
    end
    Overall_Psth_KDEfilteredSU_Action_t{nn} = cell(NT,1);
    Overall_Psth_KDEfilteredSU_Action_confb95{nn} = cell(NT,1);
    Overall_Psth_KDEfilteredSU_Action{nn} = cell(NT,1);
    for uu=1:NSU
        try
            LocalSpikes = cell2mat(SpikesTimesSU_Action{nn}(:,uu));
        catch
            LocalSpikes = cell2mat(SpikesTimesSU_Action{nn}(:,uu)');
        end
        [y,Overall_Psth_KDEfilteredSU_Action_t{nn}{uu},~, ~, ~,Overall_Psth_KDEfilteredSU_Action_confb95{nn}{uu},~]=ssvkernel(LocalSpikes,t);
        % the output of ssvkernel is a density function that sums to 1
        % multiplying by the total number of spikes gives the number of expected spike per time bin (here 10 ms)
        % dividing by the number of events gives the number of expected spike per time bin (here 10 ms) per event
        % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
        Overall_Psth_KDEfilteredSU_Action{nn}{uu} =  y * length(LocalSpikes)/Nevents * Response_samprate/1000;
    end
    
    % calculate a weighted average PSTH for each unit or tetrode across all vocalizations
    % First organize tha data into a matrix where each column represent a time
    % bin and each row a vocalization for each tetrode/unit then calculate the
    % nanmean and nanste over rows.
    Average_Psth_KDEfiltered_TAction{nn}=cell(NT,1);
    for uu=1:NT
        t=-Buffer: Bin_ms : round((max(Event_duration{nn}) + Buffer)/Bin_ms)*Bin_ms;
        Average_Psth_KDEfiltered_TAction{nn}{uu} = nan(3,length(t));
        PSTH_local = nan(Nevents,length(t));
        for ee=1:Nevents
            for tt=1:length(Psth_KDEfilteredT_Action_t{nn}{ee,uu})
                Ind = find(t==Psth_KDEfilteredT_Action_t{nn}{ee,uu}(tt));
                PSTH_local(ee,Ind) = Psth_KDEfilteredT_Action{nn}{ee,uu}(tt); %#ok<FNDSB>
            end
        end
        Average_Psth_KDEfiltered_TAction{nn}{uu}(1,:) = t;
        Average_Psth_KDEfiltered_TAction{nn}{uu}(2,:) = nanmean(PSTH_local);
        Average_Psth_KDEfiltered_TAction{nn}{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
    end
    
    Average_Psth_KDEfiltered_SUAction{nn}=cell(NSU,1);
    for uu=1:NSU
        t=-Buffer: Bin_ms : round((max(Event_duration{nn}) + Buffer)/Bin_ms)*Bin_ms;
        Average_Psth_KDEfiltered_SUAction{nn}{uu} = nan(3,length(t));
        PSTH_local = nan(Nevents,length(t));
        for vv=1:Nevents
            for tt=1:length(Psth_KDEfilteredSU_Action_t{nn}{ee,uu})
                Ind = find(t==Psth_KDEfilteredSU_Action_t{nn}{ee,uu}(tt));
                PSTH_local(vv,Ind) = Psth_KDEfilteredSU_Action{nn}{ee,uu}(tt); %#ok<FNDSB>
            end
        end
        Average_Psth_KDEfiltered_SUAction{nn}{uu}(1,:) = t;
        Average_Psth_KDEfiltered_SUAction{nn}{uu}(2,:) = nanmean(PSTH_local);
        Average_Psth_KDEfiltered_SUAction{nn}{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
    end
    
    %% Now plot Raster and KDE of the spike rate
    % Raster and average PSTH for each tetrode
    for uu=1:NT
        Fig = figure();
        subplot(2,1,1)
        for cc=1:Nevents
            hold on
            yyaxis left
            plot([0 Event_duration{nn}(cc)], cc-[0.5 0.5], '-','LineWidth',250/Nevents,'Color', [1 0.8 0.8])
            for spike=1:length(SpikesTimesT_Action{nn}{cc,uu})
                hold on
                plot(SpikesTimesT_Action{nn}{cc,uu}(spike)*ones(2,1), cc-[0.9 0.1], 'k-', 'LineWidth',1)
            end
            hold on
%             yyaxis right
%             plot(Psth_KDEfilteredT_Action_t{nn}{cc,uu}, Psth_KDEfilteredT_Action{nn}{cc,uu}/max(Psth_KDEfilteredT_Action_scalef{nn}(:,uu))+cc-1, 'r-', 'LineWidth',2)
        end
        xlabel('Time centered at production onset (ms)')
        yyaxis left
        ylim([0 Nevents+1])
        ylabel(sprintf('%s renditions', UActionText{aa}))
        title(sprintf('Raster Tetrode %d on %s for %s', uu, Date,UActionText{aa}))
        hold off
        subplot(2,1,2)
        shadedErrorBar(Average_Psth_KDEfiltered_TAction{nn}{uu}(1,:), Average_Psth_KDEfiltered_TAction{nn}{uu}(2,:), Average_Psth_KDEfiltered_TAction{nn}{uu}(3,:), {'r-', 'LineWidth',2})
        xlabel('Time centered at production onset (ms)')
        ylabel('Spike rate (/ms)')
        Ylim = get(gca, 'YLim');
        set(gca, 'YLim', [0 Ylim(2)]);
        orient(Fig,'landscape')
        Fig.PaperPositionMode = 'auto';
        set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
        set(Fig,'PaperOrientation','landscape');
        set(Fig,'PaperUnits','normalized');
        set(Fig,'PaperPosition', [0 0 1 1]);
        
        print(Fig,fullfile(Output_dir,sprintf('%s_%s_%sPSTH_Tetrode%d_%d.pdf', Date, LoggerID,UActionText{aa},uu,Buffer)),'-dpdf')
   end
    
    % KDE over all events for each tetrode
    figure()
    for uu=1:NT
        subplot(NT,1,uu)
        shadedErrorBar(Overall_Psth_KDEfilteredT_Action_t{nn}{uu}, Overall_Psth_KDEfilteredT_Action{nn}{uu}, Overall_Psth_KDEfilteredT_Action_confb95{nn}{uu}, {'k-', 'LineWidth',2})
        xlabel('Time centered at production onset (ms)')
        ylabel('Spike rate (/ms)')
        Ylim = get(gca, 'YLim');
        set(gca, 'YLim', [0 Ylim(2)]);
        title(sprintf('KDE for %s Tetrode %d on %s', UActionText{aa},uu, Date))
    end
    
    % Raster and average PSTH for each Single unit
    for uu=1:NSU
        Fig=figure();
        subplot(2,1,1)
        for cc=1:Nevents
            hold on
            yyaxis left
            plot([0 Event_duration{nn}(cc)], cc-[0.5 0.5], '-','LineWidth',250/Nevents,'Color', [1 0.8 0.8])
            for spike=1:length(SpikesTimesSU_Action{nn}{cc,uu})
                hold on
                plot(SpikesTimesSU_Action{nn}{cc,uu}(spike)*ones(2,1), cc-[0.9 0.1], 'k-', 'LineWidth',1)
            end
            hold on
%             yyaxis right
%             plot(Psth_KDEfilteredSU_Action_t{nn}{cc,uu}, Psth_KDEfilteredSU_Action{nn}{cc,uu}/max(Psth_KDEfilteredSU_Action_scalef{nn}(:,uu))+cc-1, 'r-', 'LineWidth',2)
        end
        yyaxis left
        xlabel('Time centered at production onset (ms)')
        ylim([0 Nevents+1])
        ylabel(sprintf('%s renditions', UActionText{aa}))
        title(sprintf('Raster Single Unit %d on %s for %s', uu, Date,UActionText{aa}))
        hold off
        subplot(2,1,2)
        shadedErrorBar(Average_Psth_KDEfiltered_SUAction{nn}{uu}(1,:), Average_Psth_KDEfiltered_SUAction{nn}{uu}(2,:), Average_Psth_KDEfiltered_SUAction{nn}{uu}(3,:), {'r-', 'LineWidth',2})
        xlabel('Time centered at production onset (ms)')
        ylabel('Spike rate (/ms)')
        Ylim = get(gca, 'YLim');
        set(gca, 'YLim', [0 Ylim(2)]);
        orient(Fig,'landscape')
        Fig.PaperPositionMode = 'auto';
        set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
        set(Fig,'PaperOrientation','landscape');
        set(Fig,'PaperUnits','normalized');
        set(Fig,'PaperPosition', [0 0 1 1]);
        
        print(Fig,fullfile(Output_dir,sprintf('%s_%s_%sPSTH_SU%d_%d.pdf', Date, LoggerID,UActionText{aa},uu,Buffer)),'-dpdf')
    end
    
    % KDE of rate over all events for each single unit
    figure()
    for uu=1:NSU
        subplot(NSU,1,uu)
        shadedErrorBar(Overall_Psth_KDEfilteredSU_Action_t{nn}{uu}, Overall_Psth_KDEfilteredSU_Action{nn}{uu}, Overall_Psth_KDEfilteredSU_Action_confb95{nn}{uu}, {'k-', 'LineWidth',2})
        xlabel('Time centered at production onset (ms)')
        ylabel('Spike rate (/ms)')
        Ylim = get(gca, 'YLim');
        set(gca, 'YLim', [0 Ylim(2)]);
        title(sprintf('KDE for %s Single Unit %d on %s', UActionText{aa},uu, Date))
    end
    
    
%     %Loop through events and plot LFP and RAW data for all channels for
%     %each event
%     
%     Nchannels = size(LFP_Action{nn},2);
%     for ee=1:Nevents
%         F100=figure(100);
%         for cc=1:Nchannels
%             subplot(Nchannels/2,2,cc)
%             plot(Raw_Action{nn}{ee,cc}(1:round(2*length(Raw_Action{nn}{ee,cc})/7)), '-k', 'LineWidth',2)
%             title(sprintf('Voltage trace Ch %d event %d/%d', cc-1,ee, Nevents))
%         end
%         
%         F101=figure(101);
%         for cc=1:Nchannels
%             subplot(Nchannels/2,2,cc)
%             plot(LFP_Action{nn}{ee,cc}(1:round(2*length(LFP_Action{nn}{ee,cc})/7)), '-k', 'LineWidth',2)
%             title(sprintf('LFP Ch %d event %d/%d', cc-1,ee, Nevents))
%         end
%         
%         pause(1)
%         clf(F100)
%         clf(F101)
%     end    
    
end

               
    


