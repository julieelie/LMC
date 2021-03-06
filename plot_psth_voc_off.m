function [SpikeTrains] = plot_psth_voc_off(Loggers_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, Delay, KDE_Cal,PLOT)
% AudioLoggerID = ID of the audio logger that the targeting animal is
% wearing
% NeuroLoggerID = ID of the neural logger that the targeting animal is
% wearing
% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
% (Flags(2)=1))
if nargin<8
    KDE_Cal = 0;
end
if nargin<9
    PLOT = 0;
end

% load the data
fprintf(1,'Loading data....')
load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, Delay)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'Neuro_spikes','Neuro_spikesT');
load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'FS','Piezo_wave');
% EventDir = dir(fullfile(Loggers_dir,sprintf('*%s', NeuroLoggerID(2:end)), 'extracted_data',sprintf('*%s*EVENTS.mat', Date)));
% load(fullfile(EventDir.folder, EventDir.name), 'DataDeletionOnsetOffset_usec_sync')
fprintf(1,'DONE\n')

if nargin<7
    Delay = 500;% time in ms to extract data before and after vocalization onset/offset
end


Bin_ms = 1; % size of the KDE binning
Response_samprate = 1/Bin_ms;% Sampling rate of the KDE in kHz


%% extract the spike arrival times from NeuroLoggerID for each vocalization cut of AudioLoggerID
fprintf(1, 'Gathering vocalization data...')
% centering them on the vocalization cut onset
VocInd = find(~cellfun('isempty',IndVocStartRaw_merged));
NV = length(VocInd);

% Count the number of vocalization cuts where the targeted animal is
% vocalizing vs hearing for preallocation of space
VocCall = 0;
Fns_AL = fieldnames(Piezo_wave);
if Flags(1)
    Fns_Neuro = fieldnames(Neuro_spikesT);
elseif Flags(2)
    Fns_Neuro = fieldnames(Neuro_spikes);
end
% find the vocalizations emitted by the vocalizer of interest
FocIndAudio = find(contains(Fns_AL, AudioLoggerID));
FocIndNeuro = find(contains(Fns_Neuro, NeuroLoggerID));
%     OthInd = find(~contains(Fns_AL, AudioLoggerID));
for vv=1:NV
    VocCall = VocCall + length(IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio});
end

% Now loop through calls and gather data
VocDuration = nan(1,VocCall);
if Flags(2)
    NSU = size(Neuro_spikes.(Fns_Neuro{FocIndNeuro}),2);
    SpikesTimes_VocCall = cell(VocCall,NSU);
    if KDE_Cal
        Psth_KDEfiltered_VocCall = cell(VocCall,NSU);
        Psth_KDEfiltered_VocCall_t = cell(VocCall,NSU);
        Psth_KDEfiltered_VocCall_scalef = nan(VocCall,NSU);
    end
end
if Flags(1)
    NT = size(Neuro_spikesT.(Fns_Neuro{FocIndNeuro}),2);
    SpikesTTimes_VocCall = cell(VocCall,NT);
    if KDE_Cal
        Psth_KDEfiltered_TVocCall = cell(VocCall,NT);
        Psth_KDEfiltered_TVocCall_t = cell(VocCall,NT);
        Psth_KDEfiltered_TVocCall_scalef = nan(VocCall,NT);
    end
end
% DataDeletion_VocCall = cell(VocCall,1); % This cell array contains the onset offset times (1st and 2nd column) of data deletion periods (RF Artefact) that happens during each vocalization
VocCall = 0;
Ncall = nan(NV,1);
for vv=1:NV
    Ncall(vv) = length(IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio});
    if Ncall(vv)
        for nn=1:Ncall(vv)
            VocCall = VocCall+1;
            VocDuration(VocCall) = (IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio}(nn) -IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio}(nn))/FS*1000;
            
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
            t= -round((VocDuration(VocCall) + Delay)/Bin_ms)*Bin_ms : Bin_ms : Delay;
            % Extract spikes, calculate KDE and allign everything to
            % vocalization offset
            if Flags(1)
                for uu=1:NT
                    IndT01 = logical((Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}>(IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000 - Delay)) .* (Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}<(IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000 + Delay)));
                    SpikesTTimes_VocCall{VocCall,uu} = Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}(IndT01)- IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000;
                    if KDE_Cal
                        [Psth_KDEfiltered_TVocCall{VocCall,uu}, Psth_KDEfiltered_TVocCall_t{VocCall,uu}] = kde_wrapper(SpikesTTimes_VocCall{VocCall,uu},t,Response_samprate);
                        Psth_KDEfiltered_TVocCall_scalef(VocCall,uu) = max(Psth_KDEfiltered_TVocCall{VocCall,uu});
                    end
                end
            end
            if Flags(2)
                for uu=1:NSU
                    IndSU01 = logical((Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}>(IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000 - Delay)) .* (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}<(IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000 + Delay)));
                    SpikesTimes_VocCall{VocCall,uu} = (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}(IndSU01) - IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000)';
                    if KDE_Cal
                        [Psth_KDEfiltered_VocCall{VocCall,uu}, Psth_KDEfiltered_VocCall_t{VocCall,uu}] = kde_wrapper(SpikesTimes_VocCall{VocCall,uu},t,Response_samprate);
                        Psth_KDEfiltered_VocCall_scalef(VocCall,uu) = max(Psth_KDEfiltered_VocCall{VocCall,uu});
                    end
                end
            end
        end
    end
end
fprintf(1,'DONE\n')


if sum(Ncall)
    if KDE_Cal
        fprintf(1,'Calculating average and sum KDE...')
        % calculate a weighted average PSTH for each unit or tetrode across all vocalizations
        % and single spike density estimation for all beahavioral
        % renditions
        % First organize the data into a matrix where each column represent a time
        % bin and each row a vocalization for each tetrode/unit then calculate the
        % nanmean and nanste over rows.
        if Flags(1)
%             Average_Psth_KDEfiltered_TVocCall=cell(NT,1);
            Sum_Psth_KDEfiltered_TVocCall=cell(NT,3); % This will contain the spike rate calculated from all renditions (loosing the variability between renditions) with an error estimated from the confidence interval of the spike density function
            for uu=1:NT
                t=-round((max(VocDuration) + Delay)/Bin_ms)*Bin_ms : Bin_ms : Delay;
%                 Average_Psth_KDEfiltered_TVocCall{uu} = nan(3,length(t)); % This will conain the average spike rate over renditions, it's std and the time points at which it was calculated
                PSTH_local = nan(length(VocDuration),length(t));
                for vv=1:length(VocDuration)
                    for tt=1:length(Psth_KDEfiltered_TVocCall_t{vv,uu})
                        Ind = find(t==Psth_KDEfiltered_TVocCall_t{vv,uu}(tt));
                        PSTH_local(vv,Ind) = Psth_KDEfiltered_TVocCall{vv,uu}(tt); %#ok<FNDSB>
                    end
                end
%                 Average_Psth_KDEfiltered_TVocCall{uu}(1,:) = t;
%                 Average_Psth_KDEfiltered_TVocCall{uu}(2,:) = nanmean(PSTH_local);
%                 Average_Psth_KDEfiltered_TVocCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
                
                AllSpikes_local = cell2mat(SpikesTTimes_VocCall(:,uu));
                % calculate the density estimate
                % calculate the density estimate
                [Sum_Psth_KDEfiltered_TVocCall{uu,2},Sum_Psth_KDEfiltered_TVocCall{uu,1},Sum_Psth_KDEfiltered_TVocCall{uu,3}] = kde_wrapper(AllSpikes_local,t,Response_samprate,sum(~isnan(PSTH_local)));
                fprintf(1, 'Done calculating kernel density estimate of vocalization production tetrode %d/%d\n', uu, NT);
            end
        end
        
        if Flags(2)
%             Average_Psth_KDEfiltered_VocCall=cell(NSU,1);
            Sum_Psth_KDEfiltered_VocCall=cell(NSU,3); % This will contain the spike rate calculated from all renditions (loosing the variability between renditions) with an error estimated from the confidence interval of the spike density function
            for uu=1:NSU
                t=-round((max(VocDuration) + Delay)/Bin_ms)*Bin_ms : Bin_ms : Delay; 
%                 Average_Psth_KDEfiltered_VocCall{uu} = nan(3,length(t));
                PSTH_local = nan(length(VocDuration),length(t));
                for vv=1:length(VocDuration)
                    for tt=1:length(Psth_KDEfiltered_VocCall_t{vv,uu})
                        Ind = find(t==Psth_KDEfiltered_VocCall_t{vv,uu}(tt));
                        PSTH_local(vv,Ind) = Psth_KDEfiltered_VocCall{vv,uu}(tt); %#ok<FNDSB>
                    end
                end
                
%                 Average_Psth_KDEfiltered_VocCall{uu}(1,:) = t;
%                 Average_Psth_KDEfiltered_VocCall{uu}(2,:) = nanmean(PSTH_local);
%                 Average_Psth_KDEfiltered_VocCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
                
                AllSpikes_local = cell2mat(SpikesTimes_VocCall(:,uu));
                [Sum_Psth_KDEfiltered_VocCall{uu,2},Sum_Psth_KDEfiltered_VocCall{uu,1},Sum_Psth_KDEfiltered_VocCall{uu,3}] = kde_wrapper(AllSpikes_local,t,Response_samprate,sum(~isnan(PSTH_local)));
                
                fprintf(1, 'Done calculating kernel density estimate of vocalization production SU %d/%d\n', uu, NSU);
            end
        end
    end
    
    
    
    % Now plot Raster if requested
    if PLOT
        fprintf(1,'Plotting...')
        if Flags(1)
            % We want to plot data with increasing duration of
            % vocalizations.
            [~,IDur] = sort(VocDuration);
            for uu=1:NT
                Fig=figure();
                if KDE_Cal
                    subplot(2,1,1)
                end
                
                for oo=1:VocCall
                    cc=IDur(oo);
                    hold on
                    yyaxis left
                    plot([-VocDuration(cc) 0], oo-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [1 0.8 0.8]) % vocalization
                    %         for dd=1:size(DataDeletion_VocCall{VocCall},1)
                    %             hold on
                    %             plot(DataDeletion_VocCall{VocCall}(dd,:), cc-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [0.8 0.8 0.8]) % RF artefact period
                    %         end
                    for spike=1:length(SpikesTTimes_VocCall{cc,uu})
                        hold on
                        plot(SpikesTTimes_VocCall{cc,uu}(spike)*ones(2,1), oo-[0.9 0.1], 'k-', 'LineWidth',1)
                    end
                    hold on
                    %                 yyaxis right
                    %                 plot(Psth_KDEfiltered_TVocCall_t{cc,uu}, Psth_KDEfiltered_TVocCall{cc,uu}/max(Psth_KDEfiltered_TVocCall_scalef(:,uu))+cc-1, 'r-', 'LineWidth',2)
                end
                XLIM = [-Delay-mean(VocDuration) Delay];
                xlabel('Time centered at production offset (ms)')
                yyaxis left
                ylim([0 VocCall+1])
                xlim(XLIM)
                ylabel('Vocalization production renditions')
                title(sprintf('Raster %s Tetrode %d on %s',NeuroLoggerID, uu, Date))
                hold off
                
                if KDE_Cal
                    subplot(2,1,2)
                    %                 plot(Average_Psth_KDEfiltered_TVocCall{uu}(1,:), Average_Psth_KDEfiltered_TVocCall{uu}(2,:),'r--', 'LineWidth',2)
                    %                 hold on
                    %                 plot(Sum_Psth_KDEfiltered_TVocCall{uu,1}, Sum_Psth_KDEfiltered_TVocCall{uu,2}, 'r-', 'LineWidth',2)
                    %                 hold on
                    %                 legend({'Average KDE' 'Sum KDE'}, 'AutoUpdate','off')
                    %                 shadedErrorBar(Average_Psth_KDEfiltered_TVocCall{uu}(1,:), Average_Psth_KDEfiltered_TVocCall{uu}(2,:), Average_Psth_KDEfiltered_TVocCall{uu}(3,:), {'r--', 'LineWidth',2})
                    %                 hold on
                    shadedErrorBar(Sum_Psth_KDEfiltered_TVocCall{uu,1}, Sum_Psth_KDEfiltered_TVocCall{uu,2}, Sum_Psth_KDEfiltered_TVocCall{uu,3}, {'r-', 'LineWidth',2})
                    xlim(XLIM)
                    ylim([0 max(Sum_Psth_KDEfiltered_TVocCall{uu,2})*1.3])
                    xlabel('Time centered at production offset (ms)')
                    ylabel('Spike rate (/ms)')
                    hold off
                end
                orient(Fig,'landscape')
                Fig.PaperPositionMode = 'auto';
                %             set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
                set(Fig,'PaperOrientation','landscape');
                %             set(Fig,'PaperUnits','normalized');
                %             set(Fig,'PaperPosition', [50 50 1200 800]);
                %             pause()
                if KDE_Cal
                    print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocProdPSTH_KDE_Tetrode%d_%d_off.pdf', Date, NeuroLoggerID,uu,Delay)),'-dpdf','-fillpage')
                else
                    print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocProdPSTH_Tetrode%d_%d_off.pdf', Date, NeuroLoggerID,uu,Delay)),'-dpdf','-fillpage')
                end
                close all
            end
        end
        
        if Flags(2)
            % We want to plot data with increasing duration of
            % vocalizations.
            [~,IDur] = sort(VocDuration);
            for uu=1:NSU
                Fig=figure();
                if KDE_Cal
                    subplot(2,1,1)
                end
                for oo=1:VocCall
                    cc=IDur(oo);
                    hold on
                    yyaxis left
                    plot([-VocDuration(cc) 0], oo-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [1 0.8 0.8]) % vocalizations
                    %         for dd=1:size(DataDeletion_VocCall{cc},1)
                    %             hold on
                    %             plot(DataDeletion_VocCall{cc}(dd,:), cc-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [0.8 0.8 0.8]) % RF artefact period
                    %         end
                    for spike=1:length(SpikesTimes_VocCall{cc,uu})
                        hold on
                        plot(SpikesTimes_VocCall{cc,uu}(spike)*ones(2,1), oo-[0.9 0.1], 'k-', 'LineWidth',2)
                    end
                    hold on
                    %                 yyaxis right
                    %                 plot(Psth_KDEfiltered_VocCall_t{cc,uu}, Psth_KDEfiltered_VocCall{cc,uu}/max(Psth_KDEfiltered_VocCall_scalef(:,uu))+cc-1, 'r-', 'LineWidth',2)
                end
                XLIM = [-Delay-mean(VocDuration) Delay];
                yyaxis left
                xlabel('Time centered at production offset (ms)')
                ylim([0 VocCall+1])
                xlim(XLIM)
                ylabel('Vocalization Rendition #')
                title(sprintf('Raster %s Single Unit %d on %s', NeuroLoggerID, uu, Date))
                VL = vline(0,'g:');
                VL.LineWidth=2;
                hold off
                
                if KDE_Cal
                    subplot(2,1,2)
                    %                 plot(Average_Psth_KDEfiltered_VocCall{uu}(1,:), Average_Psth_KDEfiltered_VocCall{uu}(2,:), 'r--', 'LineWidth',2)
                    %                 hold on
                    %                 plot(Sum_Psth_KDEfiltered_VocCall{uu,1}, Sum_Psth_KDEfiltered_VocCall{uu,2},'r-', 'LineWidth',2)
                    %                 legend({'Average KDE' 'Sum KDE'}, 'AutoUpdate','off')
                    %                 hold on
                    %                 shadedErrorBar(Average_Psth_KDEfiltered_VocCall{uu}(1,:), Average_Psth_KDEfiltered_VocCall{uu}(2,:), Average_Psth_KDEfiltered_VocCall{uu}(3,:), {'r--', 'LineWidth',2})
                    %                 hold on
                    shadedErrorBar(Sum_Psth_KDEfiltered_VocCall{uu,1}, Sum_Psth_KDEfiltered_VocCall{uu,2}.*10^3, Sum_Psth_KDEfiltered_VocCall{uu,3}.*10^3, {'r-', 'LineWidth',2})
                    xlim(XLIM)
                    ylim([0 max(Sum_Psth_KDEfiltered_VocCall{uu,2})*1.3*10^3])
                    xlabel('Time centered at production offset (ms)')
                    ylabel('Firing Rate (Hz)')
                    VL=vline(0,'g:');
                    VL.LineWidth=2;
                    hold off
                end
                orient(Fig,'landscape')
                Fig.PaperPositionMode = 'auto';
                %             set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
                set(Fig,'PaperOrientation','landscape');
                %             set(Fig,'PaperUnits','normalized');
                %             set(Fig,'PaperPosition', [0 0 1 1]);
                %             pause()
                if KDE_Cal
                    print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocProdPSTH_KDE_SU%d_%d_off.pdf', Date, NeuroLoggerID,uu, Delay)),'-dpdf', '-fillpage')
                else
                    print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocProdPSTH_SU%d_%d_off.pdf', Date, NeuroLoggerID,uu, Delay)),'-dpdf', '-fillpage')
                end
                close all
            end
        end
        fprintf(1, 'DONE\n')
    end
else
    fprintf('No vocalization production of the target bat on that day\n')
    Sum_Psth_KDEfiltered_TVocCall=[];
    Sum_Psth_KDEfiltered_VocCall=[];
end


%% extract the spike arrival times from NeuroLoggerID for each vocalization cut of all other audiologgers

fprintf(1, 'Gathering Hearing data...')
% centering them on the vocalization cut onset
VocInd = find(~cellfun('isempty',IndVocStartRaw_merged));
NV = length(VocInd);

% Count the number of vocalization cuts where the targeted animal is
%  hearing for preallocation of space
HearCall = 0;
% find the vocalizations heared by the vocalizer of interest
OthInd = find(~contains(Fns_AL, AudioLoggerID));
for vv=1:NV
    for ll=1:length(OthInd)
        HearCall = HearCall + length(IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)});
    end
end

% Now loop through calls and gather data
if Flags(2)
    NSU = size(Neuro_spikes.(Fns_Neuro{FocIndNeuro}),2);
    SpikesTimes_HearCall = cell(HearCall,NSU);
    if KDE_Cal
        Psth_KDEfiltered_HearCall = cell(HearCall,NSU);
        Psth_KDEfiltered_HearCall_t = cell(HearCall,NSU);
        Psth_KDEfiltered_HearCall_scalef = nan(HearCall,NSU);
    end
end
if Flags(1)
    NT = size(Neuro_spikesT.(Fns_Neuro{FocIndNeuro}),2);
    SpikesTTimes_HearCall = cell(HearCall,NT);
    if KDE_Cal
        Psth_KDEfiltered_THearCall = cell(HearCall,NT);
        Psth_KDEfiltered_THearCall_t = cell(HearCall,NT);
        Psth_KDEfiltered_THearCall_scalef = nan(HearCall,NT);
    end
end
HearDuration = nan(1,HearCall);
% DataDeletion_HearCall = cell(VocCall,1); % This cell array contains the onset offset times (1st and 2nd column) of data deletion periods (RF Artefact) that happens during each vocalization
HearOnly = nan(1,HearCall);
HearCall = 0;
for vv=1:NV
    for ll=1:length(OthInd)
        NcallHear = length(IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)});
        if NcallHear
            for nn=1:NcallHear
                HearCall = HearCall+1;
                HearDuration(HearCall) = (IndVocStopRaw_merged{VocInd(vv)}{OthInd(ll)}(nn) -IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)}(nn))/FS*1000;
                % Identify if any deletion period fall within the spike sequence for the call
%                 VocStartTransc = Voc_transc_time_refined(vv,1) + IndVocStartRaw_merged{vv}{OthInd(ll)}(nn)/FS*1000 - Delay;
%                 VocStopTransc = Voc_transc_time_refined(vv,1) + IndVocStopRaw_merged{vv}{OthInd(ll)}(nn)/FS*1000 + Delay;
%                 DelOnsetInd = find((DataDeletionOnsetOffset_usec_sync(:,1)/1000 > VocStartTransc) .* (DataDeletionOnsetOffset_usec_sync(:,1)/1000 < VocStopTransc));
%                 DelOffsetInd = find((DataDeletionOnsetOffset_usec_sync(:,2)/1000 > VocStartTransc) .* (DataDeletionOnsetOffset_usec_sync(:,2)/1000 < VocStopTransc));
%                 DelOnOffsetInd = intersect(DelOnsetInd, DelOffsetInd); % deletion periods starting and stoping within the time of the vocalization spike sequence
%                 DelOnsetInd = setdiff(DelOnsetInd, DelOnOffsetInd);  % deletion periods starting within the time of the vocalization spike sequence
%                 DelOffsetInd = setdiff(DelOffsetInd, DelOnOffsetInd);  % deletion periods stoping within the time of the vocalization spike sequence
%                 DataDeletion_HearCall{HearCall} = [DataDeletionOnsetOffset_usec_sync(DelOnOffsetInd,:)/1000; DataDeletionOnsetOffset_usec_sync(DelOnsetInd,1)/1000 VocStopTransc.*ones(length(DelOnsetInd),1) ; VocStartTransc.*ones(length(DelOffsetInd),1) DataDeletionOnsetOffset_usec_sync(DelOffsetInd,1)/1000] - VocStartTransc - Delay; % The deletion periods are time centered at vocalization onset
%                 
                t=-round((HearDuration(HearCall) + Delay)/Bin_ms)*Bin_ms : Bin_ms : Delay;
                % Is this call heard at the same time the focal bat is
                % producing a call  within the window 100ms before 100ms after a vocalization?
                OnsetSync = sum(((IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio} - Delay*FS/1000) < IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)) .* ((IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio} + Delay*FS/1000) > IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)));
                OffsetSync = sum(((IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio} - Delay*FS/1000) < IndVocStopRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)) .* ((IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio} + Delay*FS/1000) > IndVocStopRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)));
                if OnsetSync || OffsetSync % there is an overlap
                    HearOnly(HearCall) = 0;
                else % No overlap
                    HearOnly(HearCall) = 1;
                end
                if Flags(1)
                    for uu=1:NT
                        IndT01 = logical((Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}>(IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000 - Delay)) .* (Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}<(IndVocStopRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000 + Delay)));
                        SpikesTTimes_HearCall{HearCall,uu} = Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}(IndT01)- IndVocStopRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000;
                        if KDE_Cal
                            [Psth_KDEfiltered_THearCall{HearCall,uu},Psth_KDEfiltered_THearCall_t{HearCall,uu}] = kde_wrapper(SpikesTTimes_HearCall{HearCall,uu},t,Response_samprate);
                            Psth_KDEfiltered_THearCall_scalef(HearCall,uu) = max(Psth_KDEfiltered_THearCall{HearCall,uu});
                        end
                    end
                end
                if Flags(2)
                    for uu=1:NSU
                        IndSU01 = logical((Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}>(IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000 - Delay)) .* (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}<(IndVocStopRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000 + Delay)));
                        SpikesTimes_HearCall{HearCall,uu} = (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}(IndSU01) - IndVocStopRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000)';
                        if KDE_Cal
                            [Psth_KDEfiltered_HearCall{HearCall,uu},Psth_KDEfiltered_HearCall_t{HearCall,uu}] = kde_wrapper(SpikesTimes_HearCall{HearCall,uu},t,Response_samprate);
                            Psth_KDEfiltered_HearCall_scalef(HearCall,uu) = max(Psth_KDEfiltered_HearCall{HearCall,uu});
                        end
                    end
                end
            end
        end
    end
end


% calculate a weighted average PSTH for each unit or tetrode across all vocalizations
% First organize the data into a matrix where each column represent a time
% bin and each row a vocalization for each tetrode/unit then calculate the
% nanmean and nanste over rows.
% Focus only on vocalizations that were heard while the bat was not
% vocalizaing itself
if KDE_Cal
    if Flags(1) && (HearCall>1) && sum(HearOnly)
%         Average_Psth_KDEfiltered_THearCall=cell(NT,1);
        Sum_Psth_KDEfiltered_THearCall=cell(NT,3);
        HearOnlyInd = find(HearOnly);
        for uu=1:NT
            t=-round((max(HearDuration(HearOnlyInd)) + Delay)/Bin_ms)*Bin_ms: Bin_ms : Delay;
%             Average_Psth_KDEfiltered_THearCall{uu} = nan(3,length(t));
            PSTH_local = nan(length(HearOnlyInd),length(t));
            for jj=1:length(HearOnlyInd)
                vv=HearOnlyInd(jj);
                for tt=1:length(Psth_KDEfiltered_THearCall_t{vv,uu})
                    Ind = find(t==Psth_KDEfiltered_THearCall_t{vv,uu}(tt));
                    PSTH_local(vv,Ind) = Psth_KDEfiltered_THearCall{vv,uu}(tt); %#ok<FNDSB>
                end
            end
%             Average_Psth_KDEfiltered_THearCall{uu}(1,:) = t;
%             Average_Psth_KDEfiltered_THearCall{uu}(2,:) = nanmean(PSTH_local);
%             Average_Psth_KDEfiltered_THearCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
            
            AllSpikes_local = cell2mat(SpikesTTimes_HearCall(:,uu));
            % calculate the density estimate
            [Sum_Psth_KDEfiltered_THearCall{uu,2}, Sum_Psth_KDEfiltered_THearCall{uu,1}, Sum_Psth_KDEfiltered_THearCall{uu,3}] = kde_wrapper(AllSpikes_local, t,Response_samprate, (sum(~isnan(PSTH_local))));
            fprintf(1, 'Done calculating kernel density estimate of vocalization perception tetrode %d/%d\n', uu, NT);
        end
    else
        fprintf('No vocalization heard\n')
    end
    
    
    if Flags(2) && (HearCall>1) && sum(HearOnly)
%         Average_Psth_KDEfiltered_HearCall=cell(NSU,1);
        Sum_Psth_KDEfiltered_HearCall=cell(NSU,3);
        HearOnlyInd = find(HearOnly);
        for uu=1:NSU
            t=-round((max(HearDuration(HearOnlyInd)) + Delay)/Bin_ms)*Bin_ms:Bin_ms :Delay ;
%             Average_Psth_KDEfiltered_HearCall{uu} = nan(3,length(t));
            PSTH_local = nan(length(HearOnlyInd),length(t));
            for jj=1:length(HearOnlyInd)
                vv=HearOnlyInd(jj);
                for tt=1:length(Psth_KDEfiltered_HearCall_t{vv,uu})
                    Ind = find(t==Psth_KDEfiltered_HearCall_t{vv,uu}(tt));
                    PSTH_local(vv,Ind) = Psth_KDEfiltered_HearCall{vv,uu}(tt); %#ok<FNDSB>
                end
            end
%             Average_Psth_KDEfiltered_HearCall{uu}(1,:) = t;
%             Average_Psth_KDEfiltered_HearCall{uu}(2,:) = nanmean(PSTH_local);
%             Average_Psth_KDEfiltered_HearCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
            
            AllSpikes_local = cell2mat(SpikesTimes_HearCall(:,uu));
            % calculate the density estimate
            [Sum_Psth_KDEfiltered_HearCall{uu,2}, Sum_Psth_KDEfiltered_HearCall{uu,1}, Sum_Psth_KDEfiltered_HearCall{uu,3}] = kde_wrapper(AllSpikes_local,t,Response_samprate,sum(~isnan(PSTH_local)));
            
        end
    else
        fprintf('No vocalization heard\n')
    end
end
fprintf(1,'DONE\n')

% Now plot Rasters and PSTH if requested
if PLOT
    fprintf(1,'Plotting hearing data...')
    if Flags(1) && (HearCall>1) && sum(HearOnly)
        HearOnlyInd = find(HearOnly);
        HearOnlyDuration=HearDuration(HearOnlyInd);
        % We want to plot PSTH with increasing duration of vocalizations
        [~, IDurH] = sort(HearOnlyDuration);
        for uu=1:NT
            Fig=figure();
            if KDE_Cal
                subplot(2,1,1)
            end
            for hh=1:length(HearOnlyInd)
                cc = HearOnlyInd(IDurH(hh));
                yyaxis left
                hold on
                plot([-HearDuration(cc) 0], hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 1])
                %         for dd=1:size(DataDeletion_HearCall{cc},1)
                %             hold on
                %             plot(DataDeletion_HearCall{cc}(dd,:), hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 0.8]) % RF artefact period
                %         end
                for spike=1:length(SpikesTTimes_HearCall{cc,uu})
                    hold on
                    plot(SpikesTTimes_HearCall{cc,uu}(spike)*ones(2,1), hh-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
    %             yyaxis right
    %             plot(Psth_KDEfiltered_THearCall_t{cc,uu}, Psth_KDEfiltered_THearCall{cc,uu}/max(Psth_KDEfiltered_THearCall_scalef(HearOnlyInd,uu))+hh-1, 'r-', 'LineWidth',2)
            end
            XLIM = [-Delay-mean(HearDuration) Delay];
            yyaxis left
            xlabel('Time centered at hearing offset (ms)')
            ylim([0 length(HearOnlyInd)+1])
            xlim(XLIM)
            ylabel('Vocalization hearing renditions')
            title(sprintf('Raster %s Tetrode %d on %s', NeuroLoggerID, uu, Date))
            hold off
            if (HearCall>1) && KDE_Cal
                subplot(2,1,2)
    %             plot(Average_Psth_KDEfiltered_THearCall{uu}(1,:), Average_Psth_KDEfiltered_THearCall{uu}(2,:), 'b--', 'LineWidth',2)
    %             hold on
    %             plot(Sum_Psth_KDEfiltered_THearCall{uu,1}, Sum_Psth_KDEfiltered_THearCall{uu,2},'b-', 'LineWidth',2)
    %             legend({'Average KDE' 'Sum KDE'}, 'AutoUpdate', 'off')
    %             hold on
    %             shadedErrorBar(Average_Psth_KDEfiltered_THearCall{uu}(1,:), Average_Psth_KDEfiltered_THearCall{uu}(2,:), Average_Psth_KDEfiltered_THearCall{uu}(3,:), {'b--', 'LineWidth',2})
    %             hold on
                shadedErrorBar(Sum_Psth_KDEfiltered_THearCall{uu,1}, Sum_Psth_KDEfiltered_THearCall{uu,2}, Sum_Psth_KDEfiltered_THearCall{uu,3}, {'b-', 'LineWidth',2})
                xlim(XLIM)
                ylim([0 max(Sum_Psth_KDEfiltered_THearCall{uu,2})*1.3])
                xlabel('Time centered at hearing offset (ms)')
                ylabel('Spike rate (/ms)')
            end
            orient(Fig,'landscape')
            Fig.PaperPositionMode = 'auto';
    %         set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
            set(Fig,'PaperOrientation','landscape');
    %         set(Fig,'PaperUnits','normalized');
    %         set(Fig,'PaperPosition', [0 0 1 1]);
    %         pause()
            if KDE_Cal
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocHearPSTH_KDE_Tetrode%d_%d_off.pdf', Date, NeuroLoggerID,uu, Delay)),'-dpdf','-fillpage')
            else
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocHearPSTH_Tetrode%d_%d_off.pdf', Date, NeuroLoggerID,uu, Delay)),'-dpdf','-fillpage')
            end
            close all
        end
    end

    if Flags(2) && (HearCall>1) && sum(HearOnly)
        HearOnlyInd = find(HearOnly);
        HearOnlyDuration=HearDuration(HearOnlyInd);
        % We want to plot PSTH with increasing duration of vocalizations
        [~, IDurH] = sort(HearOnlyDuration);
        for uu=1:NSU
           Fig= figure();
           if KDE_Cal
               subplot(2,1,1)
           end
            for hh=1:length(HearOnlyInd)
                cc = HearOnlyInd(IDurH(hh));
                yyaxis left
                hold on
                plot([-HearDuration(cc) 0], hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 1])
                %         for dd=1:size(DataDeletion_HearCall{cc},1)
                %             hold on
                %             plot(DataDeletion_HearCall{cc}(dd,:), hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 0.8]) % RF artefact period
                %         end
                for spike=1:length(SpikesTimes_HearCall{cc,uu})
                    hold on
                    plot(SpikesTimes_HearCall{cc,uu}(spike)*ones(2,1), hh-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
    %             yyaxis right
    %             plot(Psth_KDEfiltered_HearCall_t{cc,uu}, Psth_KDEfiltered_HearCall{cc,uu}/max(Psth_KDEfiltered_HearCall_scalef(HearOnlyInd,uu))+hh-1, 'r-', 'LineWidth',2)
            end
            XLIM = [-Delay-mean(HearDuration) Delay];
            yyaxis left
            xlabel('Time centered at hearing offset (ms)')
            ylim([0 length(HearOnlyInd)+1])
            xlim(XLIM)
            ylabel('Vocalization hearing renditions')
            title(sprintf('Raster %s Single Unit %d on %s', NeuroLoggerID, uu, Date))
            hold off

            if (HearCall>1) && KDE_Cal
                subplot(2,1,2)
    %             plot(Average_Psth_KDEfiltered_HearCall{uu}(1,:), Average_Psth_KDEfiltered_HearCall{uu}(2,:), 'b-', 'LineWidth',2)
    %             hold on
    %             plot(Sum_Psth_KDEfiltered_HearCall{uu,1}, Sum_Psth_KDEfiltered_HearCall{uu,2}, Sum_Psth_KDEfiltered_HearCall{uu,3},'b-', 'LineWidth',2)
    %             legend({'Average KDE','Sum KDE'}, 'AutoUpdate', 'off')
    %             hold on
    %             shadedErrorBar(Average_Psth_KDEfiltered_HearCall{uu}(1,:), Average_Psth_KDEfiltered_HearCall{uu}(2,:), Average_Psth_KDEfiltered_HearCall{uu}(3,:), {'b-', 'LineWidth',2})
    %             hold on
                shadedErrorBar(Sum_Psth_KDEfiltered_HearCall{uu,1}, Sum_Psth_KDEfiltered_HearCall{uu,2}, Sum_Psth_KDEfiltered_HearCall{uu,3}, {'b-', 'LineWidth',2})

                xlim(XLIM)
                ylim([0 max(Sum_Psth_KDEfiltered_HearCall{uu,2})*1.3])
                xlabel('Time centered at hearing offset (ms)')
                ylabel('Spike rate (/ms)')
                hold off
            end

            orient(Fig,'landscape')
            Fig.PaperPositionMode = 'auto';
    %         set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
            set(Fig,'PaperOrientation','landscape');
    %         set(Fig,'PaperUnits','normalized');
    %         set(Fig,'PaperPosition', [0 0 1 1]);
    %         pause()
            if KDE_Cal
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocHearPSTH_KDE_SU%d_%d_off.pdf', Date, NeuroLoggerID,uu, Delay)),'-dpdf','-fillpage')
            else
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocHearPSTH_SU%d_%d_off.pdf', Date, NeuroLoggerID,uu, Delay)),'-dpdf','-fillpage')
            end
            close all
        end
    end
    fprintf(1,'DONE\n')
end

% prepare output for return
if Flags(1)
    SpikeTrains.SpikesTTimes_VocCall = SpikesTTimes_VocCall;
    SpikeTrains.SpikesTTimes_HearCall = SpikesTTimes_HearCall;
    if KDE_Cal
%         SpikeTrains.Average_Psth_KDEfiltered_TVocCall = Average_Psth_KDEfiltered_TVocCall;
        SpikeTrains.Sum_Psth_KDEfiltered_TVocCall = Sum_Psth_KDEfiltered_TVocCall;
    end
end
if Flags(2)
    SpikeTrains.SpikesTimes_VocCall = SpikesTimes_VocCall;
    SpikeTrains.SpikesTimes_HearCall = SpikesTimes_HearCall;
    if KDE_Cal
%         SpikeTrains.Average_Psth_KDEfiltered_VocCall = Average_Psth_KDEfiltered_VocCall;
        SpikeTrains.Sum_Psth_KDEfiltered_VocCall = Sum_Psth_KDEfiltered_VocCall;
    end
end
SpikeTrains.VocDuration = VocDuration;
if Flags(1)
    SpikeTrains.NT = NT;
end
if Flags(2)
    SpikeTrains.NSU = NSU;
end

if (HearCall>1) && sum(HearOnly)
    SpikeTrains.HearDuration = HearDuration;
    SpikeTrains.HearOnlyInd = HearOnlyInd;
    if Flags(1) && KDE_Cal
%         SpikeTrains.Average_Psth_KDEfiltered_THearCall = Average_Psth_KDEfiltered_THearCall;
        SpikeTrains.Sum_Psth_KDEfiltered_THearCall = Sum_Psth_KDEfiltered_THearCall;
    end
    if Flags(2) && KDE_Cal
%         SpikeTrains.Average_Psth_KDEfiltered_HearCall = Average_Psth_KDEfiltered_HearCall;
        SpikeTrains.Sum_Psth_KDEfiltered_HearCall = Sum_Psth_KDEfiltered_HearCall;
    end
else
    SpikeTrains.HearDuration = NaN;
    SpikeTrains.HearOnlyInd = NaN;
    if Flags(1) && KDE_Cal
%         SpikeTrains.Average_Psth_KDEfiltered_THearCall = NaN;
        SpikeTrains.Sum_Psth_KDEfiltered_THearCall = NaN;
    end
    if Flags(2) && KDE_Cal
%         SpikeTrains.Average_Psth_KDEfiltered_HearCall = NaN;
        SpikeTrains.Sum_Psth_KDEfiltered_HearCall = NaN;
    end
end

end
