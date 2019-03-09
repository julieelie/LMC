function [SpikeTrains] = plot_psth_behav(Loggers_dir, Date, ExpStartTime, NeuroLoggerID,BatID, Flags, MaxDur, KDE_Cal)
% AudioLoggerID = ID of the audio logger that the targeting animal is
% wearing
% NeuroLoggerID = ID of the neural logger that the targeting animal is
% wearing
% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
% (Flags(2)=1))
if nargin<8
    KDE_Cal = 0;
end

% load the data
fprintf(1,'Loading data....')
NeuroData=dir(fullfile(Loggers_dir, sprintf('%s_%s_BehavExtractData_*.mat', Date, ExpStartTime)));
load(fullfile(NeuroData.folder, NeuroData.name), 'AllActions_Time', 'AllActions_ID','UActionText','Neuro_spikes','Neuro_spikesT');
% EventDir = dir(fullfile(Loggers_dir,sprintf('*%s', NeuroLoggerID(2:end)), 'extracted_data',sprintf('*%s*EVENTS.mat', Date)));
% load(fullfile(EventDir.folder, EventDir.name), 'DataDeletionOnsetOffset_usec_sync')
fprintf(1,'DONE\n')

if nargin<7
    MaxDur = 600;% time in ms to cut behavioral actions into schrunk to constitute a PSTH
end
Response_samprate = 100;% Sampling rate of the KDE in Hz
Bin_ms = 1; % size of the KDE binning

XLIM = [0 MaxDur];
YLIM_SU = [0 0.004];
YLIM_T = [0 0.1];

%% extract the spike arrival times from NeuroLoggerID for each behavioral event and calculate a KDE at each cut if requested
if Flags(1)
    Fns_Neuro = fieldnames(Neuro_spikesT);
elseif Flags(2)
    Fns_Neuro = fieldnames(Neuro_spikes);
end
% Neuro-logger of interest
FocIndNeuro = find(contains(Fns_Neuro, NeuroLoggerID));
% Get the behavioral events
IndBehav = find(~cellfun('isempty',(Neuro_spikes.(Fns_Neuro{FocIndNeuro}))));
% Get each behavior number of cuts to allocate space
Num_Slots = nan(size(IndBehav));
for bb=1:length(IndBehav)
    Dur_local = AllActions_Time{IndBehav(bb)}(:,2)-AllActions_Time{IndBehav(bb)}(:,1);
    % actions from the focal bat
    Dur_local = Dur_local(logical(AllActions_ID{IndBehav(bb)} == BatID));
    Num_Slots(bb) = sum(floor(Dur_local/MaxDur)) + sum(mod(Dur_local, MaxDur)>0);
end


% Now loop through behaviors and gather data
if Flags(2)
    NSU = size(Neuro_spikes.(Fns_Neuro{FocIndNeuro}){IndBehav(1)},2);
    SpikesTimes_Behav = cell(length(IndBehav),1);
    if KDE_Cal
        Psth_KDEfiltered_Behav = cell(length(IndBehav),1);
        Psth_KDEfiltered_Behav_t = cell(length(IndBehav),1);
        Psth_KDEfiltered_Behav_scalef = cell(length(IndBehav),1);
    end
end
if Flags(1)
    NT = size(Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){IndBehav(1)},2);
    SpikesTTimes_Behav = cell(length(IndBehav),1);
    if KDE_Cal
        Psth_KDEfiltered_TBehav = cell(length(IndBehav),1);
        Psth_KDEfiltered_TBehav_t = cell(length(IndBehav),1);
        Psth_KDEfiltered_TBehav_scalef = cell(length(IndBehav),1);
    end
end

for bb=1:length(IndBehav)
    if Flags(2)
        SpikesTimes_Behav{bb} = cell(Num_Slots(bb),NSU);
        if KDE_Cal
            Psth_KDEfiltered_Behav{bb} = cell(Num_Slots(bb),NSU);
            Psth_KDEfiltered_Behav_t{bb} = cell(Num_Slots(bb),NSU);
            Psth_KDEfiltered_Behav_scalef{bb} = nan(Num_Slots(bb),NSU);
        end
    end
    if Flags(1)
        SpikesTTimes_Behav{bb} = cell(Num_Slots(bb),NT);
        if KDE_Cal
            Psth_KDEfiltered_TBehav{bb} = cell(Num_Slots(bb),NT);
            Psth_KDEfiltered_TBehav_t{bb} = cell(Num_Slots(bb),NT);
            Psth_KDEfiltered_TBehav_scalef{bb} = nan(Num_Slots(bb),NT);
        end
    end
    % DataDeletion_VocCall = cell(VocCall,1); % This cell array contains the onset offset times (1st and 2nd column) of data deletion periods (RF Artefact) that happens during each vocalization
    BehavCut = 0;
    FocBatInd = find(AllActions_ID{IndBehav(bb)} == BatID);
    Dur_local = AllActions_Time{IndBehav(bb)}(FocBatInd,2)-AllActions_Time{IndBehav(bb)}(FocBatInd,1);
    Nib = length(Dur_local);
    for ib=1:Nib
        NCuts = floor(Dur_local(ib)/MaxDur)+1;
        CutLimits = 0:MaxDur:(NCuts*MaxDur);
        for cc=1:NCuts
            BehavCut = BehavCut+1;
            % Identify if any deletion period fall within the spike sequence for the call
            %             VocStartTransc = Voc_transc_time_refined(vv,1) + IndVocStartRaw_merged{vv}{FocIndAudio}(nn)/FS*1000 - Delay;
            %             VocStopTransc = Voc_transc_time_refined(vv,1) + IndVocStopRaw_merged{vv}{FocIndAudio}(nn)/FS*1000 + Delay;
            %             DelOnsetInd = find((DataDeletionOnsetOffset_usec_sync(:,1)/1000 > VocStartTransc) .* (DataDeletionOnsetOffset_usec_sync(:,1)/1000 < VocStopTransc));
            %             DelOffsetInd = find((DataDeletionOnsetOffset_usec_sync(:,2)/1000 > VocStartTransc) .* (DataDeletionOnsetOffset_usec_sync(:,2)/1000 < VocStopTransc));
            %             DelOnOffsetInd = intersect(DelOnsetInd, DelOffsetInd); % deletion periods starting and stoping within the time of the vocalization spike sequence
            %             DelOnsetInd = setdiff(DelOnsetInd, DelOnOffsetInd);  % deletion periods starting within the time of the vocalization spike sequence
            %             DelOffsetInd = setdiff(DelOffsetInd, DelOnOffsetInd);  % deletion periods stoping within the time of the vocalization spike sequence
            %             DataDeletion_VocCall{VocCall} = [DataDeletionOnsetOffset_usec_sync(DelOnOffsetInd,:)/1000; DataDeletionOnsetOffset_usec_sync(DelOnsetInd,1)/1000 VocStopTransc.*ones(length(DelOnsetInd),1) ; VocStartTransc.*ones(length(DelOffsetInd),1) DataDeletionOnsetOffset_usec_sync(DelOffsetInd,1)/1000] - VocStartTransc - Delay; % The deletion periods are time centered at vocalization onset
            %
            % Calculate the t for KDE
            if cc<NCuts
                t= 0: Bin_ms : round(MaxDur/Bin_ms)*Bin_ms;
            else
                t= 0: Bin_ms : round(mod(Dur_local(ib),MaxDur)/Bin_ms)*Bin_ms;
            end
            if Flags(1)
                for uu=1:NT
                    IndT01 = logical((Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){IndBehav(bb)}{FocBatInd(ib),uu}>=CutLimits(cc)) .* (Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){IndBehav(bb)}{FocBatInd(ib),uu}<CutLimits(cc+1)));
                    SpikesTTimes_Behav{bb}{BehavCut,uu} = Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){IndBehav(bb)}{FocBatInd(ib),uu}(IndT01) - CutLimits(cc);
                    if KDE_Cal
                        if isempty(SpikesTTimes_Behav{bb}{BehavCut,uu})
                            y = ones(1, length(t))./(2*length(t));
                            Psth_KDEfiltered_TBehav{bb}{BehavCut,uu} =  y  * Response_samprate/1000;
                            Psth_KDEfiltered_TBehav_scalef{bb}(BehavCut,uu) = max(Psth_KDEfiltered_TBehav{bb}{BehavCut,uu});
                        else
                            if length(SpikesTTimes_Behav{bb}{BehavCut,uu})==1
                                y = ones(1, length(t))./(length(t));
                            else
                                % calculate the density estimate
                                [y,Psth_KDEfiltered_TBehav_t{bb}{BehavCut,uu},~]=sskernel(SpikesTTimes_Behav{bb}{BehavCut,uu},t);
                            end
                            % y is a density function that sums to 1
                            % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 10 ms)
                            % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
                            Psth_KDEfiltered_TBehav{bb}{BehavCut,uu} =  y * length(SpikesTTimes_Behav{bb}{BehavCut,uu}) * Response_samprate/1000;
                            Psth_KDEfiltered_TBehav_scalef{bb}(BehavCut,uu) = max(Psth_KDEfiltered_TBehav{bb}{BehavCut,uu});
                        end
                    end
                end
            end
            if Flags(2)
                for uu=1:NSU
                    IndSU01 = logical((Neuro_spikes.(Fns_Neuro{FocIndNeuro}){IndBehav(bb)}{FocBatInd(ib),uu}>=CutLimits(cc)) .* (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){IndBehav(bb)}{FocBatInd(ib),uu}<CutLimits(cc+1)));
                    SpikesTimes_Behav{bb}{BehavCut,uu} = Neuro_spikes.(Fns_Neuro{FocIndNeuro}){IndBehav(bb)}{FocBatInd(ib),uu}(IndSU01) - CutLimits(cc);
                    if KDE_Cal
                        if isempty(SpikesTimes_Behav{bb}{BehavCut,uu})
                            y=ones(1, length(t))./(2*length(t));
                            Psth_KDEfiltered_Behav{bb}{BehavCut,uu} =  y * Response_samprate/1000;
                            Psth_KDEfiltered_Behav_scalef{bb}(BehavCut,uu) = max(Psth_KDEfiltered_Behav{bb}{BehavCut,uu});
                        else
                            if length(SpikesTimes_Behav{bb}{BehavCut,uu})==1
                                y=ones(1,length(t))./(length(t));
                            else
                                % calculate the density estimate
                                [y,Psth_KDEfiltered_Behav_t{bb}{BehavCut,uu},~]=sskernel(SpikesTimes_Behav{bb}{BehavCut,uu},t);
                            end
                            % y is a density function that sums to 1
                            % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 10 ms)
                            % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
                            Psth_KDEfiltered_Behav{bb}{BehavCut,uu} =  y * length(SpikesTimes_Behav{bb}{BehavCut,uu}) * Response_samprate/1000;
                            Psth_KDEfiltered_Behav_scalef{bb}(BehavCut,uu) = max(Psth_KDEfiltered_Behav{bb}{BehavCut,uu});
                        end
                    end
                end
            end
        end
    end
end
fprintf(1,'DONE\n')

%% Calculating weighted PSTH
if KDE_Cal
    fprintf(1,'Calculating weighted PSTH...')
    % calculate a weighted average PSTH for each unit or tetrode across all vocalizations
    % First organize the data into a matrix where each column represent a time
    % bin and each row a vocalization for each tetrode/unit then calculate the
    % nanmean and nanste over rows.
    if Flags(1)
        Sum_Psth_KDEfiltered_TBehav=cell(NT,length(IndBehav));
        for bb=1:length(IndBehav)
            for uu=1:NT
                t=0: Bin_ms : round(MaxDur/Bin_ms)*Bin_ms;
                Sum_Psth_KDEfiltered_TBehav{uu,bb} = cell(3,1);
                NBehav = size(Psth_KDEfiltered_TBehav_t{bb},1);
                PSTH_local = nan(NBehav,length(t));
                for vv=1:NBehav
                    for tt=1:length(Psth_KDEfiltered_TBehav_t{bb}{vv,uu})
                        Ind = find(t==Psth_KDEfiltered_TBehav_t{bb}{vv,uu}(tt));
                        PSTH_local(vv,Ind) = Psth_KDEfiltered_TBehav{bb}{vv,uu}(tt); %#ok<FNDSB>
                    end
                end
                
                AllSpikes_local = cell2mat(SpikesTTimes_Behav{bb}(:,uu));
                % calculate the density estimate
                if isempty(AllSpikes_local)
                    y=ones(1,length(t))./(2*length(t));
                    Sum_Psth_KDEfiltered_TBehav{uu,bb}{1} = t;
                    Sum_Psth_KDEfiltered_TBehav{uu,bb}{2} =  y ./(sum(~isnan(PSTH_local)))* Response_samprate/1000;
                    Sum_Psth_KDEfiltered_TBehav{uu,bb}{3} = nan(2, length(t));
                else
                    if length(AllSpikes_local)==1
                        y=ones(1,length(t))./length(t);
                        Sum_Psth_KDEfiltered_TBehav{uu,bb}{1} = t;
                    else
                        % calculate the density estimate
                        [y,Sum_Psth_KDEfiltered_TBehav{uu,bb}{1},OptW,~,~,bconf95]=sskernel(AllSpikes_local,t);
                    end
                    fprintf(1, 'Done calculating kernel density estimate for %s events tetrode %d/%d, using kernel bandwidth= %f\n', UActionText{IndBehav(bb)}, uu, NT, OptW);
                    % y is a density function that sums to 1
                    % Multiplying by the total number of spikes gives the number of expecting spike per time bin for all behavioral events (here 10 ms)
                    % Dividing by the number of behavioral events per time bin
                    % gives the number of expected spikes per behavioral event
                    % Multiplying by the response sampling rate in kHz gives the expected spike rate to one behavioral event in spike/ms
                    Sum_Psth_KDEfiltered_TBehav{uu,bb}{2} =  y * length(AllSpikes_local) ./(sum(~isnan(PSTH_local))) * Response_samprate/1000;
                    Sum_Psth_KDEfiltered_TBehav{uu,bb}{3} =  abs(flipud(bconf95) .* length(AllSpikes_local) ./(sum(~isnan(PSTH_local))) .* Response_samprate ./1000 - repmat(Sum_Psth_KDEfiltered_TBehav{uu,bb}{2},2,1));
                end
           end
        end
    end
    
    if Flags(2)
        Sum_Psth_KDEfiltered_Behav=cell(NSU,length(IndBehav));
        for bb=1:length(IndBehav)
            for uu=1:NSU
                t=0: Bin_ms : round(MaxDur/Bin_ms)*Bin_ms;
                Sum_Psth_KDEfiltered_Behav{uu,bb} = cell(3,1);
                NBehav = size(Psth_KDEfiltered_Behav_t{bb},1);
                PSTH_local = nan(NBehav,length(t));
                for vv=1:NBehav
                    for tt=1:length(Psth_KDEfiltered_Behav_t{bb}{vv,uu})
                        Ind = find(t==Psth_KDEfiltered_Behav_t{bb}{vv,uu}(tt));
                        PSTH_local(vv,Ind) = Psth_KDEfiltered_Behav{bb}{vv,uu}(tt); %#ok<FNDSB>
                    end
                end
                
                AllSpikes_local = cell2mat(SpikesTimes_Behav{bb}(:,uu));
                % calculate the density estimate
                if isempty(AllSpikes_local)
                    y=ones(1,length(t))./(2*length(t));
                    Sum_Psth_KDEfiltered_Behav{uu,bb}{1} = t;
                    Sum_Psth_KDEfiltered_Behav{uu,bb}{2} =  y ./(sum(~isnan(PSTH_local)))* Response_samprate/1000;
                    Sum_Psth_KDEfiltered_Behav{uu,bb}{3} = nan(2, length(t));
                else
                    if length(AllSpikes_local)==1
                        y=ones(1,length(t))./length(t);
                        Sum_Psth_KDEfiltered_Behav{uu,bb}{1} = t;
                    else
                        % calculate the density estimate
                        [y,Sum_Psth_KDEfiltered_Behav{uu,bb}{1},OptW,~,~,bconf95]=sskernel(AllSpikes_local,t);
                    end
                    fprintf(1, 'Done calculating kernel density estimate for %s events tetrode %d/%d, using kernel bandwidth= %f\n', UActionText{IndBehav(bb)}, uu, NT, OptW);
                    % y is a density function that sums to 1
                    % Multiplying by the total number of spikes gives the number of expecting spike per time bin for all behavioral events (here 10 ms)
                    % Dividing by the number of behavioral events per time bin
                    % gives the number of expected spikes per behavioral event
                    % Multiplying by the response sampling rate in kHz gives the expected spike rate to one behavioral event in spike/ms
                    Sum_Psth_KDEfiltered_Behav{uu,bb}{2} =  y * length(AllSpikes_local) ./(sum(~isnan(PSTH_local))) * Response_samprate/1000;
                    Sum_Psth_KDEfiltered_Behav{uu,bb}{3} =  abs(fliup(bconf95) .* length(AllSpikes_local) ./(sum(~isnan(PSTH_local))) .* Response_samprate ./1000 - repmat(Sum_Psth_KDEfiltered_Behav{uu,bb}{2},2,1));
                end
            end
        end
    end
    fprintf(1,'DONE\n')
end

%% Now plotting PSTH
fprintf(1, 'Plotting PSTH\n')
% Now plot Raster
if Flags(1)
    for bb=1:length(IndBehav)
        NBehav = size(SpikesTTimes_Behav{bb},1);
        for uu=1:NT
            Fig=figure();
            if KDE_Cal
                subplot(2,1,1)
            end
            
            for cc=1:NBehav
                hold on
%                 yyaxis left
                for spike=1:length(SpikesTTimes_Behav{bb}{cc,uu})
                    hold on
                    plot(SpikesTTimes_Behav{bb}{cc,uu}(spike)*ones(2,1), cc-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
                %                 yyaxis right
                %                 plot(Psth_KDEfiltered_TVocCall_t{cc,uu}, Psth_KDEfiltered_TVocCall{cc,uu}/max(Psth_KDEfiltered_TVocCall_scalef(:,uu))+cc-1, 'r-', 'LineWidth',2)
            end
            xlabel('Time (continuous during behavior, in ms)')
%             yyaxis left
            ylim([0 NBehav+1])
            xlim(XLIM)
            ylabel('Renditions')
            title(sprintf('%s Raster %s Tetrode %d on %s', UActionText{IndBehav(bb)},NeuroLoggerID, uu, Date))
            hold off
            
            if KDE_Cal
                subplot(2,1,2)
                shadedErrorBar(Sum_Psth_KDEfiltered_TBehav{uu,bb}{1}, Sum_Psth_KDEfiltered_TBehav{uu,bb}{2}, flipud(Sum_Psth_KDEfiltered_TBehav{uu,bb}{3})-Sum_Psth_KDEfiltered_TBehav{uu,bb}{2}, {'r-', 'LineWidth',2})
                xlim(XLIM)
                ylim(YLIM_T)
                xlabel('Time (ms)')
                ylabel('Spike rate (/ms)')
            end
            orient(Fig,'landscape')
            Fig.PaperPositionMode = 'auto';
            %             set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
            set(Fig,'PaperOrientation','landscape');
            %             set(Fig,'PaperUnits','normalized');
            %             set(Fig,'PaperPosition', [50 50 1200 800]);
            %             pause()
            if KDE_Cal
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_%sPSTH_KDE_Tetrode%d.pdf', Date, NeuroLoggerID,UActionText{IndBehav(bb)},uu)),'-dpdf','-fillpage')
            else
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_%sPSTH_Tetrode%d.pdf', Date, NeuroLoggerID,UActionText{IndBehav(bb)},uu)),'-dpdf','-fillpage')
            end
            close all
        end
    end
end

if Flags(2)
    for bb=1:length(IndBehav)
        NBehav = size(SpikesTimes_Behav{bb},1);
        for uu=1:NSU
            Fig=figure();
            if KDE_Cal
                subplot(2,1,1)
            end
            for cc=1:NBehav
                hold on
                for spike=1:length(SpikesTimes_Behav{bb}{cc,uu})
                    hold on
                    plot(SpikesTimes_Behav{bb}{cc,uu}(spike)*ones(2,1), cc-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
                %                 yyaxis right
                %                 plot(Psth_KDEfiltered_VocCall_t{cc,uu}, Psth_KDEfiltered_VocCall{cc,uu}/max(Psth_KDEfiltered_VocCall_scalef(:,uu))+cc-1, 'r-', 'LineWidth',2)
            end
            yyaxis left
            xlabel('Time (continuous during behavior, in ms)')
            ylim([0 NBehav+1])
            xlim(XLIM)
            ylabel('Renditions')
            title(sprintf('%s Raster %s Single Unit %d on %s', UActionText{IndBehav(bb)},NeuroLoggerID, uu, Date))
            hold off
            
            if KDE_Cal
                subplot(2,1,2)
                shadedErrorBar(Sum_Psth_KDEfiltered_Behav{uu,bb}{1}, Sum_Psth_KDEfiltered_Behav{uu,bb}{2}, flipud(Sum_Psth_KDEfiltered_Behav{uu,bb}{3})-Sum_Psth_KDEfiltered_Behav{uu,bb}{2}, {'r-', 'LineWidth',2})
                xlim(XLIM)
                ylim(YLIM_SU)
                xlabel('Time (ms)')
                ylabel('Spike rate (/ms)')
            end
            orient(Fig,'landscape')
            Fig.PaperPositionMode = 'auto';
            %             set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
            set(Fig,'PaperOrientation','landscape');
            %             set(Fig,'PaperUnits','normalized');
            %             set(Fig,'PaperPosition', [0 0 1 1]);
            %             pause()
            if KDE_Cal
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_%sPSTH_KDE_SU%d.pdf', Date, NeuroLoggerID,UActionText{IndBehav(bb)},uu)),'-dpdf', '-fillpage')
            else
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_%sPSTH_SU%d.pdf', Date, NeuroLoggerID,UActionText{IndBehav(bb)},uu)),'-dpdf', '-fillpage')
            end
            close all
        end
    end
end
fprintf(1, 'DONE\n')
if Flags(2)
    SpikeTrains.SpikesTimes_Behav = SpikesTimes_Behav;
    if KDE_Cal
        SpikeTrains.Sum_Psth_KDEfiltered_TBehav = Sum_Psth_KDEfiltered_TBehav;
    end
end
if Flags(1)
    SpikeTrains.SpikesTTimes_Behav = SpikesTTimes_Behav;
    if KDE_Cal
        SpikeTrains.Sum_Psth_KDEfiltered_Behav = Sum_Psth_KDEfiltered_Behav;
    end
end
SpikeTrains.UActionBehav = UActionText(IndBehav);
end
