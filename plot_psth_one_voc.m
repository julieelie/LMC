function plot_psth_one_voc(SpikeTimes, Raw_wave, FS)

% FigureSize = [1 1 30 20];

Response_samprate = 100;% Sampling rate of the KDE in Hz
Bin_ms = 1; % size of the KDE binning

DB_noise = 60; % Noise threshold for the spectrogram colormap
FHigh_spec = 90000; % Max frequency (Hz) for the raw data spectrogram
BandPassFilter = [1000 90000];% Frequency bands chosen for digital signal processing
Fhigh_power =20; % Frequency upper bound for calculating the envelope (time running RMS)
Fs_env = 1000; % Sample frequency of the enveloppe

%% plot the spike arrival times for each sound section of Raw_wave
NV = length(Raw_wave);

% Now loop through calls and gather data
VocDuration = nan(1,NV); % Duration of each sound extract in ms
NU = size(SpikeTimes,2);
Psth_KDEfiltered_VocCall = cell(NV,NU);
Psth_KDEfiltered_VocCall_t = cell(NV,NU);
Psth_KDEfiltered_VocCall_scalef = nan(NV,NU);


for vv=1:NV
    VocDuration(vv) = length(Raw_wave{vv}/FS)*1000; %#ok<IDISVAR,USENS>
    
    % Calculate the t for KDE
    t=0: Bin_ms : round(VocDuration(vv));
    
    for uu=1:NU
        % calculate the density estimate
        [y,Psth_KDEfiltered_VocCall_t{vv,uu},~]=ssvkernel(SpikeTimes{vv,uu},t);
        % y is a density function that sums to 1
        % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 10 ms)
        % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
        Psth_KDEfiltered_VocCall{vv,uu} =  y * length(SpikeTimes{vv,uu}) * Response_samprate/1000;
        Psth_KDEfiltered_VocCall_scalef(vv,uu) = max(Psth_KDEfiltered_VocCall{vv,uu});
    end
end

XLIM = [0 max(VocDuration)];

% calculate a weighted average PSTH for each unit or tetrode across all vocalizations
% First organize tha data into a matrix where each column represent a time
% bin and each row a vocalization for each tetrode/unit then calculate the
% nanmean and nanste over rows.
Average_Psth_KDEfiltered_VocCall=cell(NU,1);
for uu=1:NU
    t=0: Bin_ms : round(max(VocDuration)/Bin_ms)*Bin_ms;
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


% design filters of raw ambient recording, bandpass and low pass which was
% used for the cross correlation
[z,p,k] = butter(6,BandPassFilter/(FS/2),'bandpass');
sos_raw_band = zp2sos(z,p,k);

% Now plot Raster
for vv=1:NV
    % plot the spectrogram of the raw_wave
    % bandpass filter the ambient mic recording
    Filt_RawVoc = filtfilt(sos_raw_band,1,Raw_wave{vv});
    Amp_env_Mic = running_rms(Filt_RawVoc, FS, Fhigh_power, Fs_env);
    Fig=figure();
    subplot(2,1,1)
    [~] = spec_only_bats(Filt_RawVoc, FS, DB_noise, FHigh_spec);
    hold on
    yyaxis right
    plot((1:length(Amp_env_Mic{vv}))/Fs_env*1000, Amp_env_Mic{vv}, 'r-', 'LineWidth',2)
    ylabel('Amplitude')
    title(sprintf('Ambient Microphone Voc %d/%d',vv,NV))
    
    for uu=1:NU
        subplot(2,1,2)
        hold on
        yyaxis left
        for spike=1:length(SpikeTimes{uu,vv})
            hold on
            plot(SpikeTimes{uu,vv}(spike)*ones(2,1), uu-[0.9 0.1], 'k-', 'LineWidth',1)
        end
        hold on
        yyaxis right
        plot(Psth_KDEfiltered_VocCall_t{uu,vv}, Psth_KDEfiltered_VocCall{uu,vv}/max(Psth_KDEfiltered_VocCall_scalef(:,vv))+uu-1, 'r-', 'LineWidth',2)
    end
    xlabel('Time centered at productionsection onset (ms)')
    yyaxis left
    ylim([0 NU+1])
    xlim(XLIM)
    ylabel('Units')
    title(sprintf('Raster sound %d on %s', vv, Date))
    hold off
    
    pause()
    
%     orient(Fig,'landscape')
%     Fig.PaperPositionMode = 'auto';
%     set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
%     set(Fig,'PaperOrientation','landscape');
%     set(Fig,'PaperUnits','normalized');
%     set(Fig,'PaperPosition', [0 0 1 1]);
%     print(Fig,fullfile(Loggers_dir,sprintf('%s_Voc%d_PSTH_AllUnits.pdf', Date,vv)),'-dpdf')
end
%             figure()
%             shadedErrorBar(Average_Psth_KDEfiltered_VocCall{uu}(1,:), Average_Psth_KDEfiltered_VocCall{uu}(2,:), Average_Psth_KDEfiltered_VocCall{uu}(3,:), {'r-', 'LineWidth',2})
%             xlim(XLIM)
%             ylim(YLIM_T)
%             xlabel('Time centered at sound onset (ms)')
%             ylabel('Spike rate (/ms)')

end
