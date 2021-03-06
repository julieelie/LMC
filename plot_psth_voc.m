function [SpikeTrains] = plot_psth_voc(Loggers_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, Delay, KDE_Cal,PLOT,RthreshInd)
% AudioLoggerID = ID of the audio logger that the targeting animal is
% wearing
% NeuroLoggerID = ID of the neural logger that the targeting animal is
% wearing
% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
% (Flags(2)=1))
% All KDE are calculated in spike/ms or kHz
if nargin<8
    KDE_Cal = 0;
end
if nargin<9
    PLOT = 0;
end
if nargin<10
    RthreshInd = [];
end

% load the data
fprintf(1,'Loading data....')
DataFile = fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, Delay));
load(DataFile, 'IndVocStartRaw_merged', 'IndVocStopRaw_merged','Neuro_spikes_Baseline','BSL_transc_time_refined');
if Flags(1)
    load(DataFile, 'Neuro_spikesT');
    if ~isempty(RthreshInd)
        load(DataFile,'Neuro_spikesTDeNoiseInd')
    end
end
if Flags(2)
    load(DataFile, 'Neuro_spikes');
end

load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'FS','Piezo_wave','Raw_wave');
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
VocWave = cell(1,VocCall);
Voc_BSLDuration = nan(1,VocCall);
if Flags(2)
    NSU = size(Neuro_spikes.(Fns_Neuro{FocIndNeuro}),2);
    SpikesTimes_VocBaseline = cell(VocCall,NSU);
    SpikesTimes_VocCall = cell(VocCall,NSU);
    if KDE_Cal
        Psth_KDEfiltered_VocCall = cell(VocCall,NSU);
        Psth_KDEfiltered_VocCall_t = cell(VocCall,NSU);
        Psth_KDEfiltered_VocCall_scalef = nan(VocCall,NSU);
        Psth_KDEfiltered_VocBaseline = cell(VocCall,NSU);
        Psth_KDEfiltered_VocBaseline_t = cell(VocCall,NSU);
        Psth_KDEfiltered_VocBaseline_scalef = nan(VocCall,NSU);
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
            IndOn = IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio}(nn) - Delay*FS/1000;
            IndOff = IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio}(nn)+Delay*FS/1000;
            DurWave_local = length(Raw_wave{VocInd(vv)});
            WL = Raw_wave{VocInd(vv)}( max(1,IndOn):min(DurWave_local, IndOff));
            if IndOn<1
                WL = [zeros(abs(IndOn),1); WL];
            end
            if IndOff>DurWave_local
                WL = [WL ; zeros(IndOff-DurWave_local,1)];
            end
            VocWave{VocCall} = WL; % contains the vocalization with the same delay before/after as the neural response
            Voc_BSLDuration(VocCall) = BSL_transc_time_refined(VocInd(vv),2)-BSL_transc_time_refined(VocInd(vv),1);
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
            t=-Delay: Bin_ms : round((VocDuration(VocCall) + Delay)/Bin_ms)*Bin_ms;
            t_baseline=0: Bin_ms : Voc_BSLDuration(VocCall);
            if Flags(1)
                for uu=1:NT
                    IndT01 = logical((Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}>(IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000 - Delay)) .* (Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}<(IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000 + Delay)));
                    if ~isempty(RthreshInd)
                        IndT01 = logical(IndT01 .* Neuro_spikesTDeNoiseInd.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}{RthreshInd});
                    end
                    SpikesTTimes_VocCall{VocCall,uu} = Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}(IndT01)- IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000;
                    if KDE_Cal
                        [Psth_KDEfiltered_TVocCall{VocCall,uu}, Psth_KDEfiltered_TVocCall_t{VocCall,uu}] = kde_wrapper(SpikesTTimes_VocCall{VocCall,uu},t,Response_samprate);
                        Psth_KDEfiltered_TVocCall_scalef(VocCall,uu) = max(Psth_KDEfiltered_TVocCall{VocCall,uu});
                    end
                end
            end
            if Flags(2)
                for uu=1:NSU
                    IndSU01 = logical((Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}>(IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000 - Delay)) .* (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}<(IndVocStopRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000 + Delay)));
                    SpikesTimes_VocCall{VocCall,uu} = (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}(IndSU01) - IndVocStartRaw_merged{VocInd(vv)}{FocIndAudio}(nn)/FS*1000)';
                    if KDE_Cal
                        [Psth_KDEfiltered_VocCall{VocCall,uu}, Psth_KDEfiltered_VocCall_t{VocCall,uu}] = kde_wrapper(SpikesTimes_VocCall{VocCall,uu},t,Response_samprate);
                        Psth_KDEfiltered_VocCall_scalef(VocCall,uu) = max(Psth_KDEfiltered_VocCall{VocCall,uu});
                        
                        % calculate a KDE for the baseline period
                        if ~isnan(Voc_BSLDuration(VocCall)) % if isnan: No baseline period could be taken before the onset of the vocalization
                            SpikesTimes_VocBaseline{VocCall,uu} = Neuro_spikes_Baseline.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu};
                            [Psth_KDEfiltered_VocBaseline{VocCall,uu}, Psth_KDEfiltered_VocBaseline_t{VocCall,uu}] = kde_wrapper(SpikesTimes_VocBaseline{VocCall,uu},t_baseline,Response_samprate);
                            Psth_KDEfiltered_VocBaseline_scalef(VocCall,uu) = max(Psth_KDEfiltered_VocBaseline{VocCall,uu});
                        end
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
                t=-Delay: Bin_ms : round((max(VocDuration) + Delay)/Bin_ms)*Bin_ms;
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
                [Sum_Psth_KDEfiltered_TVocCall{uu,2},Sum_Psth_KDEfiltered_TVocCall{uu,1},Sum_Psth_KDEfiltered_TVocCall{uu,3}] = kde_wrapper(AllSpikes_local,t,Response_samprate,sum(~isnan(PSTH_local)));
                fprintf(1, 'Done calculating kernel density estimate of vocalization production tetrode %d/%d\n', uu, NT);
            end
        end
        
        if Flags(2)
%             Average_Psth_KDEfiltered_VocCall=cell(NSU,1);
            Sum_Psth_KDEfiltered_VocCall=cell(NSU,3); % This will contain the spike rate calculated from all renditions (loosing the variability between renditions) with an error estimated from the confidence interval of the spike density function
            Sum_Psth_KDEfiltered_VocBaseline=cell(NSU,3); % This will contain the spike rate calculated from all baseline (loosing the variability between renditions) with an error estimated from the confidence interval of the spike density function
            for uu=1:NSU
                t=-Delay: Bin_ms : round((max(VocDuration) + Delay)/Bin_ms)*Bin_ms;
                t_baseline = 0:Bin_ms:max(Voc_BSLDuration);
%                 Average_Psth_KDEfiltered_VocCall{uu} = nan(3,length(t));
                PSTH_local = nan(length(VocDuration),length(t));
                for vv=1:length(VocDuration)
                    for tt=1:length(Psth_KDEfiltered_VocCall_t{vv,uu})
                        Ind = find(t==Psth_KDEfiltered_VocCall_t{vv,uu}(tt));
                        PSTH_local(vv,Ind) = Psth_KDEfiltered_VocCall{vv,uu}(tt); %#ok<FNDSB>
                    end
                end
                Baseline_data = nan(length(Voc_BSLDuration),length(t_baseline));
                for vv=1:length(Voc_BSLDuration)
                    for tt=1:length(Psth_KDEfiltered_VocBaseline_t{vv,uu})
                        Ind = find(t_baseline==Psth_KDEfiltered_VocBaseline_t{vv,uu}(tt));
                        Baseline_data(vv,Ind) = Psth_KDEfiltered_VocBaseline{vv,uu}(tt); %#ok<FNDSB>
                    end
                end
%                 Average_Psth_KDEfiltered_VocCall{uu}(1,:) = t;
%                 Average_Psth_KDEfiltered_VocCall{uu}(2,:) = nanmean(PSTH_local);
%                 Average_Psth_KDEfiltered_VocCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
                
                AllSpikes_local = cell2mat(SpikesTimes_VocCall(:,uu));
                [Sum_Psth_KDEfiltered_VocCall{uu,2},Sum_Psth_KDEfiltered_VocCall{uu,1},Sum_Psth_KDEfiltered_VocCall{uu,3}] = kde_wrapper(AllSpikes_local,t,Response_samprate,sum(~isnan(PSTH_local)));
                

                AllSpikesBaseline_local = cell2mat(SpikesTimes_VocBaseline(:,uu)');
                USpikeT = unique(AllSpikesBaseline_local);
                CountSpikeT = histc(AllSpikesBaseline_local,USpikeT);
                RepSpike = USpikeT(CountSpikeT>1);
                for rs=1:length(RepSpike)
                    Ind_sr = find(AllSpikesBaseline_local == RepSpike(rs));
                    AllSpikesBaseline_local(Ind_sr) = AllSpikesBaseline_local(Ind_sr) + [0 ones(1,length(Ind_sr)-1)]; % add 1ms to the spike that are duplicated so the calculation of the PSTH can go through and is still balanced compared to the PSTh during vocalizations
                end
                [Sum_Psth_KDEfiltered_VocBaseline{uu,2},Sum_Psth_KDEfiltered_VocBaseline{uu,1},Sum_Psth_KDEfiltered_VocBaseline{uu,3}] = kde_wrapper(AllSpikesBaseline_local,t_baseline,Response_samprate,sum(~isnan(Baseline_data)));
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
                    plot([0 VocDuration(cc)], oo-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [1 0.8 0.8]) % vocalization
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
                XLIM = [-Delay mean(VocDuration)+Delay];
                xlabel('Time centered at production onset (ms)')
%                 yyaxis left
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
                    xlabel('Time centered at production onset (ms)')
                    ylabel('Spike rate (/ms or kHz)')
                    hold off
                end
                orient(Fig,'landscape')
                Fig.PaperPositionMode = 'auto';
                %             set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
                set(Fig,'PaperOrientation','landscape');
                %             set(Fig,'PaperUnits','normalized');
                %             set(Fig,'PaperPosition', [50 50 1200 800]);
                %             pause()
                if KDE_Cal && ~isempty(RthreshInd)
                    Filename_local = sprintf('%s_%s_VocProdPSTH_KDE_Tetrode%d_%d_DeNoise%d.pdf', Date, NeuroLoggerID,uu, Delay,RthreshInd);
                elseif KDE_Cal && isempty(RthreshInd)
                    Filename_local = sprintf('%s_%s_VocProdPSTH_KDE_Tetrode%d_%d.pdf', Date, NeuroLoggerID,uu, Delay);
                else
                    Filename_local = sprintf('%s_%s_VocProdPSTH_Tetrode%d_%d.pdf', Date, NeuroLoggerID,uu, Delay);
                end
                print(Fig,fullfile(Loggers_dir,Filename_local),'-dpdf','-fillpage')
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
%                     yyaxis left
                    plot([0 VocDuration(cc)], oo-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [1 0.8 0.8]) % vocalizations
                    %         for dd=1:size(DataDeletion_VocCall{cc},1)
                    %             hold on
                    %             plot(DataDeletion_VocCall{cc}(dd,:), cc-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [0.8 0.8 0.8]) % RF artefact period
                    %         end
                    for spike=1:length(SpikesTimes_VocCall{cc,uu})
                        hold on
%                         plot(SpikesTimes_VocCall{cc,uu}(spike)*ones(2,1), oo-[0.9 0.1], 'k-', 'LineWidth',2)
                        plot(SpikesTimes_VocCall{cc,uu}(spike), oo, 'ko','MarkerSize',3,'MarkerFaceColor', 'k')
                    end
                    hold on
                    %                 yyaxis right
                    %                 plot(Psth_KDEfiltered_VocCall_t{cc,uu}, Psth_KDEfiltered_VocCall{cc,uu}/max(Psth_KDEfiltered_VocCall_scalef(:,uu))+cc-1, 'r-', 'LineWidth',2)
                end
                XLIM = [-Delay mean(VocDuration)+Delay];
%                 yyaxis left
                xlabel('Time centered at production onset (ms)')
                ylim([0 VocCall+1])
                set(gca,'YTick', [0 50 100 150], 'YTickLabel', [0 50 100 150], 'XTick',[],'XTickLabel',[])
                xlim(XLIM)
                ylabel('Vocalization Rendition #')
                title(sprintf('Raster %s Single Unit %d on %s', NeuroLoggerID, uu, Date))
                VL = vline(0,'b:');
                VL.LineWidth=2;
                hold off
                
                if KDE_Cal
%                     figure()
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
                    ylim([0 max(Sum_Psth_KDEfiltered_VocCall{uu,2}(logical((Sum_Psth_KDEfiltered_VocCall{uu,1}>XLIM(1)).*(Sum_Psth_KDEfiltered_VocCall{uu,1}<XLIM(2))))).*1.3*10^3])
                    xlabel('Time (ms)')
                    ylabel('Firing Rate (Hz)')
%                     set(gca,'YTick',[0 5 10 15],'YTickLabel',[0 5 10 15])
                    VL=vline(0,'b:');
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
                    Filename_local = sprintf('%s_%s_VocProdPSTH_KDE_SU%d_%d.pdf', Date, NeuroLoggerID,uu, Delay);
                else
                    Filename_local = sprintf('%s_%s_VocProdPSTH_SU%d_%d.pdf', Date, NeuroLoggerID,uu, Delay);
                end
                print(Fig,fullfile(Loggers_dir,Filename_local),'-dpdf', '-fillpage')
                close all
            end
        end
        fprintf(1, 'DONE\n')
    end
else
    fprintf('No vocalization production of the target bat on that day\n')
    Sum_Psth_KDEfiltered_TVocCall=[];
    Sum_Psth_KDEfiltered_VocCall=[];
    Sum_Psth_KDEfiltered_VocBaseline=[];
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
    SpikesTimes_HearBaseline = cell(HearCall,NSU);
    if KDE_Cal
        Psth_KDEfiltered_HearCall = cell(HearCall,NSU);
        Psth_KDEfiltered_HearCall_t = cell(HearCall,NSU);
        Psth_KDEfiltered_HearCall_scalef = nan(HearCall,NSU);
        Psth_KDEfiltered_HearBaseline = cell(HearCall,NSU);
        Psth_KDEfiltered_HearBaseline_t = cell(HearCall,NSU);
        Psth_KDEfiltered_HearBaseline_scalef = nan(HearCall,NSU);
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
Hear_BSLDuration = nan(1,HearCall);
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
                Hear_BSLDuration(HearCall) = BSL_transc_time_refined(VocInd(vv),2)-BSL_transc_time_refined(VocInd(vv),1);
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
                t=-Delay: Bin_ms : round((HearDuration(HearCall) + Delay)/Bin_ms)*Bin_ms;
                t_baseline=0: Bin_ms : round((Hear_BSLDuration(HearCall))/Bin_ms)*Bin_ms;
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
                        if ~isempty(RthreshInd)
                            IndT01 = logical(IndT01 .* Neuro_spikesTDeNoiseInd.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}{RthreshInd});
                        end
                        SpikesTTimes_HearCall{HearCall,uu} = Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}(IndT01)- IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000;
                        if KDE_Cal
                            [Psth_KDEfiltered_THearCall{HearCall,uu},Psth_KDEfiltered_THearCall_t{HearCall,uu}] = kde_wrapper(SpikesTTimes_HearCall{HearCall,uu},t,Response_samprate);
                            Psth_KDEfiltered_THearCall_scalef(HearCall,uu) = max(Psth_KDEfiltered_THearCall{HearCall,uu});
                        end
                    end
                end
                if Flags(2)
                    for uu=1:NSU
                        IndSU01 = logical((Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}>(IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000 - Delay)) .* (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}<(IndVocStopRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000 + Delay)));
                        SpikesTimes_HearCall{HearCall,uu} = (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu}(IndSU01) - IndVocStartRaw_merged{VocInd(vv)}{OthInd(ll)}(nn)/FS*1000)';
                        if KDE_Cal
                            [Psth_KDEfiltered_HearCall{HearCall,uu},Psth_KDEfiltered_HearCall_t{HearCall,uu}] = kde_wrapper(SpikesTimes_HearCall{HearCall,uu},t,Response_samprate);
                            Psth_KDEfiltered_HearCall_scalef(HearCall,uu) = max(Psth_KDEfiltered_HearCall{HearCall,uu});
                                
                            SpikesTimes_HearBaseline{HearCall,uu} = Neuro_spikes_Baseline.(Fns_Neuro{FocIndNeuro}){VocInd(vv),uu};
                            [Psth_KDEfiltered_HearBaseline{HearCall,uu},Psth_KDEfiltered_HearBaseline_t{HearCall,uu}] = kde_wrapper(SpikesTimes_HearBaseline{HearCall,uu},t_baseline,Response_samprate);
                            Psth_KDEfiltered_HearBaseline_scalef(HearCall,uu) = max(Psth_KDEfiltered_HearBaseline{HearCall,uu});
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
            t=-Delay: Bin_ms : round((max(HearDuration(HearOnlyInd)) + Delay)/Bin_ms)*Bin_ms;
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
        Sum_Psth_KDEfiltered_HearBaseline = cell(NSU,3);
        HearOnlyInd = find(HearOnly);
        for uu=1:NSU
            t=-Delay: Bin_ms : round((max(HearDuration(HearOnlyInd)) + Delay)/Bin_ms)*Bin_ms;
            t_baseline=0: Bin_ms : round((max(Hear_BSLDuration(HearOnlyInd)))/Bin_ms)*Bin_ms;
%             Average_Psth_KDEfiltered_HearCall{uu} = nan(3,length(t));
            PSTH_local = nan(length(HearOnlyInd),length(t));
            for jj=1:length(HearOnlyInd)
                vv=HearOnlyInd(jj);
                for tt=1:length(Psth_KDEfiltered_HearCall_t{vv,uu})
                    Ind = find(t==Psth_KDEfiltered_HearCall_t{vv,uu}(tt));
                    PSTH_local(vv,Ind) = Psth_KDEfiltered_HearCall{vv,uu}(tt); %#ok<FNDSB>
                end
            end
            Baseline_data = nan(length(HearOnlyInd),length(t_baseline));
            for jj=1:length(HearOnlyInd)
                vv=HearOnlyInd(jj);
                for tt=1:length(Psth_KDEfiltered_HearBaseline_t{vv,uu})
                    Ind = find(t_baseline==Psth_KDEfiltered_HearBaseline_t{vv,uu}(tt));
                    Baseline_data(vv,Ind) = Psth_KDEfiltered_HearBaseline{vv,uu}(tt); %#ok<FNDSB>
                end
            end
%             Average_Psth_KDEfiltered_HearCall{uu}(1,:) = t;
%             Average_Psth_KDEfiltered_HearCall{uu}(2,:) = nanmean(PSTH_local);
%             Average_Psth_KDEfiltered_HearCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
            
             AllSpikes_local = cell2mat(SpikesTimes_HearCall(:,uu));
            % calculate the density estimate
            [Sum_Psth_KDEfiltered_HearCall{uu,2}, Sum_Psth_KDEfiltered_HearCall{uu,1}, Sum_Psth_KDEfiltered_HearCall{uu,3}] = kde_wrapper(AllSpikes_local,t,Response_samprate,sum(~isnan(PSTH_local)));
            
            AllSpikesBaseline_local = cell2mat(SpikesTimes_HearBaseline(:,uu)');
            USpikeT = unique(AllSpikesBaseline_local);
            CountSpikeT = histc(AllSpikesBaseline_local,USpikeT);
            RepSpike = USpikeT(CountSpikeT>1);
            for rs=1:length(RepSpike)
                Ind_sr = find(AllSpikesBaseline_local == RepSpike(rs));
                AllSpikesBaseline_local(Ind_sr) = AllSpikesBaseline_local(Ind_sr) + [0 ones(1,length(Ind_sr)-1)]; % add 1ms to the spike that are duplicated so the calculation of the PSTH can go through and is still balanced compared to the PSTh during vocalizations
            end
            [Sum_Psth_KDEfiltered_HearBaseline{uu,2}, Sum_Psth_KDEfiltered_HearBaseline{uu,1}, Sum_Psth_KDEfiltered_HearBaseline{uu,3}] = kde_wrapper(AllSpikesBaseline_local,t_baseline,Response_samprate,sum(~isnan(Baseline_data)));
            
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
                plot([0 HearDuration(cc)], hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 1])
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
            XLIM = [-Delay mean(HearDuration)+Delay];
            yyaxis left
            xlabel('Time centered at hearing onset (ms)')
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
                xlabel('Time centered at hearing onset (ms)')
                ylabel('Spike rate (/ms) or kHz')
            end
            orient(Fig,'landscape')
            Fig.PaperPositionMode = 'auto';
    %         set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
            set(Fig,'PaperOrientation','landscape');
    %         set(Fig,'PaperUnits','normalized');
    %         set(Fig,'PaperPosition', [0 0 1 1]);
    %         pause()
            if KDE_Cal && ~isempty(RthreshInd)
                Filename_local = sprintf('%s_%s_VocHearPSTH_KDE_Tetrode%d_%d_Denoise%d.pdf', Date, NeuroLoggerID,uu, Delay,RthreshInd);
            elseif KDE_Cal && isempty(RthreshInd)
                Filename_local = sprintf('%s_%s_VocHearPSTH_KDE_Tetrode%d_%d.pdf', Date, NeuroLoggerID,uu, Delay);
            else
                Filename_local = sprintf('%s_%s_VocHearPSTH_Tetrode%d_%d.pdf', Date, NeuroLoggerID,uu, Delay);
            end
            print(Fig,fullfile(Loggers_dir,Filename_local),'-dpdf','-fillpage')
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
                plot([0 HearDuration(cc)], hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 1])
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
            XLIM = [-Delay mean(HearDuration)+Delay];
            yyaxis left
            xlabel('Time centered at hearing onset (ms)')
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
                xlabel('Time centered at hearing onset (ms)')
                ylabel('Spike rate (/ms or kHz)')
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
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocHearPSTH_KDE_SU%d_%d.pdf', Date, NeuroLoggerID,uu, Delay)),'-dpdf','-fillpage')
            else
                print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_VocHearPSTH_SU%d_%d.pdf', Date, NeuroLoggerID,uu, Delay)),'-dpdf','-fillpage')
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
    SpikeTrains.NT = NT;
end
if Flags(2)
    SpikeTrains.SpikesTimes_VocCall = SpikesTimes_VocCall;
    SpikeTrains.SpikesTimes_HearCall = SpikesTimes_HearCall;
    SpikeTrains.SpikesTimes_HearBaseline = SpikesTimes_HearBaseline;
    SpikeTrains.SpikesTimes_VocBaseline = SpikesTimes_VocBaseline;
    if KDE_Cal
%         SpikeTrains.Average_Psth_KDEfiltered_VocCall = Average_Psth_KDEfiltered_VocCall;
        SpikeTrains.Sum_Psth_KDEfiltered_VocCall = Sum_Psth_KDEfiltered_VocCall;
        SpikeTrains.Sum_Psth_KDEfiltered_VocBaseline = Sum_Psth_KDEfiltered_VocBaseline;
    end
    SpikeTrains.NSU = NSU;
end
SpikeTrains.VocDuration = VocDuration;
SpikeTrains.VocWave =  VocWave;


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
        SpikeTrains.Sum_Psth_KDEfiltered_HearBaseline = Sum_Psth_KDEfiltered_HearBaseline;
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
        SpikeTrains.Sum_Psth_KDEfiltered_HearBaseline = NaN;
    end
end

end
