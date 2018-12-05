function [AllActions, UActionText]= plot_psth_event(ExtData_dir, Date,AllActions, UActionText, Buffer, LoggerID, Output_dir)
% AllActions is a cell array of length the number of actions containing
% each 2 columns giving the onset and offset times in ms of behavioral
% events, constructed by get_logger_data_event.m
if nargin<5
    Buffer = 500; %ms
end

AllSep = strfind(ExtData_dir, filesep);
if nargin<6
    if length(ExtData_dir)==AllSep(end)
        LoggerID = ExtData_dir((AllSep(end-2)+1):(AllSep(end-1)-1));
    else
        LoggerID = ExtData_dir((AllSep(end-1)+1):(AllSep(end)-1));
    end
end
if nargin<7
    if length(ExtData_dir)==AllSep(end)
        Output_dir = ExtData_dir(1:(AllSep(end-2)));
    else
        Output_dir = ExtData_dir(1:(AllSep(end-1)));
    end
end
fprintf(1, 'Figures and data are saved to %s\n', Output_dir)
MaxEventDur = 1000; %ms

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
        
        print(Fig,fullfile(Output_dir,sprintf('%s_%s_%sPSTH_Tetrode%d_%d.pdf', Date(3:end), LoggerID,UActionText{aa},uu,Buffer)),'-dpdf')
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
        
        print(Fig,fullfile(Output_dir,sprintf('%s_%s_%sPSTH_SU%d_%d.pdf', Date(3:end), LoggerID,UActionText{aa},uu,Buffer)),'-dpdf')
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

               
    


