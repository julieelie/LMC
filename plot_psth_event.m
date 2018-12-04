function plot_psth_event(Loggers_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, Delay)
% AudioLoggerID = ID of the audio logger that the targeting animal is
% wearing
% NeuroLoggerID = ID of the neural logger that the targeting animal is
% wearing
% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
% (Flags(2)=1))



% load the data
load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged','Neuro_LFP','Neuro_spikes','Neuro_spikesT', 'FS','Piezo_FS','Piezo_wave','VocFilename','Voc_transc_time_refined');
EventDir = dir(fullfile(Loggers_dir,sprintf('l%s', NeuroLoggerID(2:end)), 'extracted_data',sprintf('*%s*EVENTS.mat', Date)));
load(fullfile(EventDir.folder, EventDir.name), 'DataDeletionOnsetOffset_usec_sync')
if nargin<7
    Delay = 500;% time in ms to extract data before and after vocalization onset/offset
end
Response_samprate = 100;% Sampling rate of the KDE in Hz
Bin_ms = 10; % size of the KDE binning

XLIM = [-Delay 1000 + Delay];
YLIM_SU = [0 0.002];
YLIM_T = [0 0.005];

%% extract the spike arrival times from NeuroLoggerID for each vocalization cut of AudioLoggerID
% centering them on the vocalization cut onset
NV = length(VocFilename);

% Count the number of vocalization cuts where the targeted animal is
% vocalizing vs hearing for preallocation of space
VocCall = 0;
Fns_AL = fieldnames(Piezo_wave);
Fns_Neuro = fieldnames(Neuro_spikes);
 % find the vocalizations emitted by the vocalizer of interest
FocIndAudio = find(contains(Fns_AL, AudioLoggerID));
FocIndNeuro = find(contains(Fns_Neuro, NeuroLoggerID));
%     OthInd = find(~contains(Fns_AL, AudioLoggerID));
for vv=1:NV
   VocCall = VocCall + length(IndVocStartRaw_merged{vv}{FocIndAudio});%#ok<IDISVAR,USENS>
end

% Now loop through calls and gather data
VocDuration = nan(1,VocCall);
if Flags(2)
    NSU = size(Neuro_spikes.(Fns_Neuro{FocIndNeuro}),2);
    SpikesTimes_VocCall = cell(VocCall,NSU);
    Psth_KDEfiltered_VocCall = cell(VocCall,NSU);
    Psth_KDEfiltered_VocCall_t = cell(VocCall,NSU);
    Psth_KDEfiltered_VocCall_scalef = nan(VocCall,NSU);
end
if Flags(1)
    NT = size(Neuro_spikesT.(Fns_Neuro{FocIndNeuro}),2);
    SpikesTTimes_VocCall = cell(VocCall,NT);
    Psth_KDEfiltered_TVocCall = cell(VocCall,NT);
    Psth_KDEfiltered_TVocCall_t = cell(VocCall,NT);
    Psth_KDEfiltered_TVocCall_scalef = nan(VocCall,NT);
end
% DataDeletion_VocCall = cell(VocCall,1); % This cell array contains the onset offset times (1st and 2nd column) of data deletion periods (RF Artefact) that happens during each vocalization
VocCall = 0;
Ncall = nan(NV,1);
for vv=1:NV
    Ncall(vv) = length(IndVocStartRaw_merged{vv}{FocIndAudio});
    if Ncall(vv)
        for nn=1:Ncall(vv)
            VocCall = VocCall+1;
            VocDuration(VocCall) = (IndVocStopRaw_merged{vv}{FocIndAudio}(nn) -IndVocStartRaw_merged{vv}{FocIndAudio}(nn))/FS*1000; %#ok<IDISVAR,USENS>
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
            if Flags(1)
                for uu=1:NT
                    IndT01 = logical((Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){vv,uu}>(IndVocStartRaw_merged{vv}{FocIndAudio}(nn)/FS*1000 - Delay)) .* (Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){vv,uu}<(IndVocStopRaw_merged{vv}{FocIndAudio}(nn)/FS*1000 + Delay)));
                    SpikesTTimes_VocCall{VocCall,uu} = Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){vv,uu}(IndT01)- IndVocStartRaw_merged{vv}{FocIndAudio}(nn)/FS*1000;
                    % calculate the density estimate
                    [y,Psth_KDEfiltered_TVocCall_t{VocCall,uu},~]=ssvkernel(SpikesTTimes_VocCall{VocCall,uu},t);
                    % y is a density function that sums to 1
                    % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 10 ms)
                    % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
                    Psth_KDEfiltered_TVocCall{VocCall,uu} =  y * length(SpikesTTimes_VocCall{VocCall,uu}) * Response_samprate/1000;
                    Psth_KDEfiltered_TVocCall_scalef(VocCall,uu) = max(Psth_KDEfiltered_TVocCall{VocCall,uu});
                end
            end
            if Flags(2)
                for uu=1:NSU
                    IndSU01 = logical((Neuro_spikes.(Fns_Neuro{FocIndNeuro}){vv,uu}>(IndVocStartRaw_merged{vv}{FocIndAudio}(nn)/FS*1000 - Delay)) .* (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){vv,uu}<(IndVocStopRaw_merged{vv}{FocIndAudio}(nn)/FS*1000 + Delay)));
                    SpikesTimes_VocCall{VocCall,uu} = Neuro_spikes.(Fns_Neuro{FocIndNeuro}){vv,uu}(IndSU01) - IndVocStartRaw_merged{vv}{FocIndAudio}(nn)/FS*1000;
                    % calculate the density estimate
                    [y,Psth_KDEfiltered_VocCall_t{VocCall,uu},~]=ssvkernel(SpikesTimes_VocCall{VocCall,uu},t);
                    % y is a density function that sums to 1
                    % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 10 ms)
                    % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
                    Psth_KDEfiltered_VocCall{VocCall,uu} =  y * length(SpikesTimes_VocCall{VocCall,uu}) * Response_samprate/1000;
                    Psth_KDEfiltered_VocCall_scalef(VocCall,uu) = max(Psth_KDEfiltered_VocCall{VocCall,uu});
                end
            end
        end
    end
end


if sum(Ncall)
        % calculate a weighted average PSTH for each unit or tetrode across all vocalizations
    % First organize tha data into a matrix where each column represent a time
    % bin and each row a vocalization for each tetrode/unit then calculate the
    % nanmean and nanste over rows.
    if Flags(1)
        Average_Psth_KDEfiltered_TVocCall=cell(NT,1);
        for uu=1:NT
            t=-Delay: Bin_ms : round((max(VocDuration) + Delay)/Bin_ms)*Bin_ms;
            Average_Psth_KDEfiltered_TVocCall{uu} = nan(3,length(t));
            PSTH_local = nan(length(VocDuration),length(t));
            for vv=1:length(VocDuration)
                for tt=1:length(Psth_KDEfiltered_TVocCall_t{vv,uu})
                    Ind = find(t==Psth_KDEfiltered_TVocCall_t{vv,uu}(tt));
                    PSTH_local(vv,Ind) = Psth_KDEfiltered_TVocCall{vv,uu}(tt); %#ok<FNDSB>
                end
            end
            Average_Psth_KDEfiltered_TVocCall{uu}(1,:) = t;
            Average_Psth_KDEfiltered_TVocCall{uu}(2,:) = nanmean(PSTH_local);
            Average_Psth_KDEfiltered_TVocCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
        end
    end

    if Flags(2)
        Average_Psth_KDEfiltered_VocCall=cell(NSU,1);
        for uu=1:NSU
            t=-Delay: Bin_ms : round((max(VocDuration) + Delay)/Bin_ms)*Bin_ms;
            Average_Psth_KDEfiltered_VocCall{uu} = nan(3,length(t));
            PSTH_local = nan(length(VocDuration),length(t));
            for vv=1:length(VocDuration)
                for tt=1:length(Psth_KDEfiltered_VocCall_t{vv,uu})
                    Ind = find(t==Psth_KDEfiltered_VocCall_t{vv,uu}(tt));
                    PSTH_local(vv,Ind) = Psth_KDEfiltered_VocCall{vv,uu}(tt); %#ok<FNDSB>
                end
            end
            Average_Psth_KDEfiltered_VocCall{uu}(1,:) = t;
            Average_Psth_KDEfiltered_VocCall{uu}(2,:) = nanmean(PSTH_local);
            Average_Psth_KDEfiltered_VocCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
        end
    end



    % Now plot Raster
    if Flags(1)
        for uu=1:NT
            figure()
            subplot(2,1,1)
            for cc=1:VocCall
                hold on
                yyaxis left
                plot([0 VocDuration(cc)], cc-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [1 0.8 0.8]) % vocalization
                %         for dd=1:size(DataDeletion_VocCall{VocCall},1)
                %             hold on
                %             plot(DataDeletion_VocCall{VocCall}(dd,:), cc-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [0.8 0.8 0.8]) % RF artefact period
                %         end
                for spike=1:length(SpikesTTimes_VocCall{cc,uu})
                    hold on
                    plot(SpikesTTimes_VocCall{cc,uu}(spike)*ones(2,1), cc-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
%                 yyaxis right
%                 plot(Psth_KDEfiltered_TVocCall_t{cc,uu}, Psth_KDEfiltered_TVocCall{cc,uu}/max(Psth_KDEfiltered_TVocCall_scalef(:,uu))+cc-1, 'r-', 'LineWidth',2)
            end
            xlabel('Time centered at production onset (ms)')
            yyaxis left
            ylim([0 VocCall+1])
            xlim(XLIM)
            ylabel('Vocalization production renditions')
            title(sprintf('Raster Tetrode %d on %s', uu, Date))
            hold off

            subplot(2,1,2)
            shadedErrorBar(Average_Psth_KDEfiltered_TVocCall{uu}(1,:), Average_Psth_KDEfiltered_TVocCall{uu}(2,:), Average_Psth_KDEfiltered_TVocCall{uu}(3,:), {'r-', 'LineWidth',2})
            xlim(XLIM)
            ylim(YLIM_T)
            xlabel('Time centered at production onset (ms)')
            ylabel('Spike rate (/ms)')
            set(gcf, 'Position', get(0, 'Screensize'));
            saveas(gcf,fullfile(Loggers_dir,sprintf('%s_%s_VocProdPSTH_Tetrode%d_%d.pdf', Date, NeuroLoggerID,uu,Delay)),'pdf')
        end
    end

    if Flags(2)
        for uu=1:NSU
            figure()
            subplot(2,1,1)
            for cc=1:VocCall
                hold on
                yyaxis left
                plot([0 VocDuration(cc)], cc-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [1 0.8 0.8]) % vocalizations
                %         for dd=1:size(DataDeletion_VocCall{cc},1)
                %             hold on
                %             plot(DataDeletion_VocCall{cc}(dd,:), cc-[0.5 0.5], '-','LineWidth',250/VocCall,'Color', [0.8 0.8 0.8]) % RF artefact period
                %         end
                for spike=1:length(SpikesTimes_VocCall{cc,uu})
                    hold on
                    plot(SpikesTimes_VocCall{cc,uu}(spike)*ones(2,1), cc-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
%                 yyaxis right
%                 plot(Psth_KDEfiltered_VocCall_t{cc,uu}, Psth_KDEfiltered_VocCall{cc,uu}/max(Psth_KDEfiltered_VocCall_scalef(:,uu))+cc-1, 'r-', 'LineWidth',2)
            end
            yyaxis left
            xlabel('Time centered at production onset (ms)')
            ylim([0 VocCall+1])
            xlim(XLIM)
            ylabel('Vocalization production renditions')
            title(sprintf('Raster Single Unit %d on %s', uu, Date))
            hold off

            subplot(2,1,2)
            shadedErrorBar(Average_Psth_KDEfiltered_VocCall{uu}(1,:), Average_Psth_KDEfiltered_VocCall{uu}(2,:), Average_Psth_KDEfiltered_VocCall{uu}(3,:), {'r-', 'LineWidth',2})
            xlim(XLIM)
            ylim(YLIM_SU)
            xlabel('Time centered at production onset (ms)')
            ylabel('Spike rate (/ms)')
            set(gcf, 'Position', get(0, 'Screensize'));
            saveas(gcf,fullfile(Loggers_dir,sprintf('%s_%s_VocProdPSTH_SU%d_%d.pdf', Date, NeuroLoggerID,uu, Delay)),'pdf')
        end
    end
else
    fprintf('No vocalization production of the target bat on that day\n')
end

%% extract the spike arrival times from NeuroLoggerID for each vocalization cut of all other audiologgers
% centering them on the vocalization cut onset
NV = length(VocFilename);

% Count the number of vocalization cuts where the targeted animal is
%  hearing for preallocation of space
HearCall = 0;
% find the vocalizations heared by the vocalizer of interest
OthInd = find(~contains(Fns_AL, AudioLoggerID));
for vv=1:NV
    for ll=1:length(OthInd)
        HearCall = HearCall + length(IndVocStartRaw_merged{vv}{OthInd(ll)});
    end
end

% Now loop through calls and gather data
if Flags(2)
    NSU = size(Neuro_spikes.(Fns_Neuro{FocIndNeuro}),2);
    SpikesTimes_HearCall = cell(HearCall,NSU);
    Psth_KDEfiltered_HearCall = cell(HearCall,NSU);
    Psth_KDEfiltered_HearCall_t = cell(HearCall,NSU);
    Psth_KDEfiltered_HearCall_scalef = nan(HearCall,NSU);
end
if Flags(1)
    NT = size(Neuro_spikesT.(Fns_Neuro{FocIndNeuro}),2);
    SpikesTTimes_HearCall = cell(HearCall,NT);
    Psth_KDEfiltered_THearCall = cell(HearCall,NT);
    Psth_KDEfiltered_THearCall_t = cell(HearCall,NT);
    Psth_KDEfiltered_THearCall_scalef = nan(HearCall,NT);
end
HearDuration = nan(1,HearCall);
% DataDeletion_HearCall = cell(VocCall,1); % This cell array contains the onset offset times (1st and 2nd column) of data deletion periods (RF Artefact) that happens during each vocalization
HearOnly = nan(1,HearCall);
HearCall = 0;
for vv=1:NV
    for ll=1:length(OthInd)
        NcallHear = length(IndVocStartRaw_merged{vv}{OthInd(ll)});
        if NcallHear
            for nn=1:NcallHear
                HearCall = HearCall+1;
                HearDuration(HearCall) = (IndVocStopRaw_merged{vv}{OthInd(ll)}(nn) -IndVocStartRaw_merged{vv}{OthInd(ll)}(nn))/FS*1000;
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
                % Is this call heard at the same time the focal bat is
                % producing a call  within the window 100ms before 100ms after a vocalization?
                OnsetSync = sum(((IndVocStartRaw_merged{vv}{FocIndAudio} - Delay*FS/1000) < IndVocStartRaw_merged{vv}{OthInd(ll)}(nn)) .* ((IndVocStopRaw_merged{vv}{FocIndAudio} + Delay*FS/1000) > IndVocStartRaw_merged{vv}{OthInd(ll)}(nn)));
                OffsetSync = sum(((IndVocStartRaw_merged{vv}{FocIndAudio} - Delay*FS/1000) < IndVocStopRaw_merged{vv}{OthInd(ll)}(nn)) .* ((IndVocStopRaw_merged{vv}{FocIndAudio} + Delay*FS/1000) > IndVocStopRaw_merged{vv}{OthInd(ll)}(nn)));
                if OnsetSync || OffsetSync % there is an overlap
                    HearOnly(HearCall) = 0;
                else % No overlap
                    HearOnly(HearCall) = 1;
                end
                if Flags(1)
                    for uu=1:NT
                        IndT01 = logical((Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){vv,uu}>(IndVocStartRaw_merged{vv}{OthInd(ll)}(nn)/FS*1000 - Delay)) .* (Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){vv,uu}<(IndVocStopRaw_merged{vv}{OthInd(ll)}(nn)/FS*1000 + Delay)));
                        SpikesTTimes_HearCall{HearCall,uu} = Neuro_spikesT.(Fns_Neuro{FocIndNeuro}){vv,uu}(IndT01)- IndVocStartRaw_merged{vv}{OthInd(ll)}(nn)/FS*1000;
                        % calculate the density estimate
                        [y,Psth_KDEfiltered_THearCall_t{HearCall,uu},~]=ssvkernel(SpikesTTimes_HearCall{HearCall,uu},t);
                        % y is a density function that sums to 1
                        % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 10 ms)
                        % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
                        Psth_KDEfiltered_THearCall{HearCall,uu} =  y * length(SpikesTTimes_HearCall{HearCall,uu}) * Response_samprate/1000;
                        Psth_KDEfiltered_THearCall_scalef(HearCall,uu) = max(Psth_KDEfiltered_THearCall{HearCall,uu});
                    end
                end
                if Flags(2)
                    for uu=1:NSU
                        IndSU01 = logical((Neuro_spikes.(Fns_Neuro{FocIndNeuro}){vv,uu}>(IndVocStartRaw_merged{vv}{OthInd(ll)}(nn)/FS*1000 - Delay)) .* (Neuro_spikes.(Fns_Neuro{FocIndNeuro}){vv,uu}<(IndVocStopRaw_merged{vv}{OthInd(ll)}(nn)/FS*1000 + Delay)));
                        SpikesTimes_HearCall{HearCall,uu} = Neuro_spikes.(Fns_Neuro{FocIndNeuro}){vv,uu}(IndSU01) - IndVocStartRaw_merged{vv}{OthInd(ll)}(nn)/FS*1000;
                        % calculate the density estimate
                        [y,Psth_KDEfiltered_HearCall_t{HearCall,uu},~]=ssvkernel(SpikesTimes_HearCall{HearCall,uu},t);
                        % y is a density function that sums to 1
                        % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 10 ms)
                        % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
                        Psth_KDEfiltered_HearCall{HearCall,uu} =  y * length(SpikesTimes_HearCall{HearCall,uu}) * Response_samprate/1000;
                        Psth_KDEfiltered_HearCall_scalef(HearCall,uu) = max(Psth_KDEfiltered_HearCall{HearCall,uu});
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
if Flags(1)
    Average_Psth_KDEfiltered_THearCall=cell(NT,1);
    HearOnlyInd = find(HearOnly);
    for uu=1:NT
        t=-Delay: Bin_ms : round((max(HearDuration(HearOnlyInd)) + Delay)/Bin_ms)*Bin_ms;
        Average_Psth_KDEfiltered_THearCall{uu} = nan(3,length(t));
        PSTH_local = nan(length(HearOnlyInd),length(t));
        for jj=1:length(HearOnlyInd)
            vv=HearOnlyInd(jj);
            for tt=1:length(Psth_KDEfiltered_THearCall_t{vv,uu})
                Ind = find(t==Psth_KDEfiltered_THearCall_t{vv,uu}(tt));
                PSTH_local(vv,Ind) = Psth_KDEfiltered_THearCall{vv,uu}(tt); %#ok<FNDSB>
            end
        end
        Average_Psth_KDEfiltered_THearCall{uu}(1,:) = t;
        Average_Psth_KDEfiltered_THearCall{uu}(2,:) = nanmean(PSTH_local);
        Average_Psth_KDEfiltered_THearCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
    end
end

if Flags(2)
    Average_Psth_KDEfiltered_HearCall=cell(NSU,1);
    HearOnlyInd = find(HearOnly);
    for uu=1:NSU
        t=-Delay: Bin_ms : round((max(HearDuration(HearOnlyInd)) + Delay)/Bin_ms)*Bin_ms;
        Average_Psth_KDEfiltered_HearCall{uu} = nan(3,length(t));
        PSTH_local = nan(length(HearOnlyInd),length(t));
        for jj=1:length(HearOnlyInd)
            vv=HearOnlyInd(jj);
            for tt=1:length(Psth_KDEfiltered_HearCall_t{vv,uu})
                Ind = find(t==Psth_KDEfiltered_HearCall_t{vv,uu}(tt));
                PSTH_local(vv,Ind) = Psth_KDEfiltered_HearCall{vv,uu}(tt); %#ok<FNDSB>
            end
        end
        Average_Psth_KDEfiltered_HearCall{uu}(1,:) = t;
        Average_Psth_KDEfiltered_HearCall{uu}(2,:) = nanmean(PSTH_local);
        Average_Psth_KDEfiltered_HearCall{uu}(3,:) = nanstd(PSTH_local)./(sum(~isnan(PSTH_local))).^0.5;
    end
end

% Now plot Rasters and PSTH
if Flags(1)
    for uu=1:NT
        figure()
        subplot(2,1,1)
        for hh=1:length(HearOnlyInd)
            cc = HearOnlyInd(hh);
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
        yyaxis left
        xlabel('Time centered at hearing onset (ms)')
        ylim([0 length(HearOnlyInd)+1])
        xlim(XLIM)
        ylabel('Vocalization hearing renditions')
        title(sprintf('Raster Tetrode %d on %s', uu, Date))
        hold off
        
        subplot(2,1,2)
        shadedErrorBar(Average_Psth_KDEfiltered_THearCall{uu}(1,:), Average_Psth_KDEfiltered_THearCall{uu}(2,:), Average_Psth_KDEfiltered_THearCall{uu}(3,:), {'b-', 'LineWidth',2})
        xlim(XLIM)
        ylim(YLIM_T)
        xlabel('Time centered at hearing onset (ms)')
        ylabel('Spike rate (/ms)')
        
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,fullfile(Loggers_dir,sprintf('%s_%s_VocHearPSTH_Tetrode%d_%d.pdf', Date, NeuroLoggerID,uu, Delay)),'pdf')
    end
end

if Flags(2)
    for uu=1:NSU
        figure()
        subplot(2,1,1)
        for hh=1:length(HearOnlyInd)
            cc = HearOnlyInd(hh);
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
        yyaxis left
        xlabel('Time centered at hearing onset (ms)')
        ylim([0 length(HearOnlyInd)+1])
        xlim(XLIM)
        ylabel('Vocalization hearing renditions')
        title(sprintf('Raster Single Unit %d on %s', uu, Date))
        hold off
        
        subplot(2,1,2)
        shadedErrorBar(Average_Psth_KDEfiltered_HearCall{uu}(1,:), Average_Psth_KDEfiltered_HearCall{uu}(2,:), Average_Psth_KDEfiltered_HearCall{uu}(3,:), {'b-', 'LineWidth',2})
        xlim(XLIM)
        ylim(YLIM_SU)
        xlabel('Time centered at hearing onset (ms)')
        ylabel('Spike rate (/ms)')
        
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,fullfile(Loggers_dir,sprintf('%s_%s_VocHearPSTH_SU%d_%d.pdf', Date, NeuroLoggerID,uu, Delay)),'pdf')
    end
end
