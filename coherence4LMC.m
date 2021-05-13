function Coherence = coherence4LMC(Cell, Self, Session, FeatureName, Trill, Nvoc,BootstrapType, Average, Delay)
% Self is a boolean 0: vocalizations from other bats (auditory coherency) 1: vocalizations from self: motor coherency
% Session = 'Operant' or 'Free'
% FeatureName = 'amp' or 'SpectralMean' or 'sal' The feature on which the
% coherece should be calculated
% Trill=0; % Set to 1 to only do calculations on Trill calls, set to 0 to do calculations on all calls
% here we calculate the coherency between the neural response and the
% acoustic features
if nargin<4
    FeatureName = 'amp';
end
if nargin<5
    Trill=0;
end
if nargin<7
    BootstrapType=0;
end

if nargin<6
    Nvoc=0;
end

if nargin<8
    Average = 0; %1 calculate X as the weighted average over all input X
end

if nargin<9
    Delay=200; % The segment of data taken into account is -Delay ms before the vocalization onset and +200ms after the vocalization offset
end

PlotCoherenceFig = 1; % To plot the result of coherence calculation for each cell
StimXYDataPlot=0;
TR=2; % 2ms is chosen as the Time resolution for the neural data
Fs = 1/(TR*10^-3); % the data is then sampled at the optimal frequency given the neural time resolution choosen
% find the closest power of 2 for the number of FFT window points that
% correspond to the Nyquist limit
Nyquist = Fs * 0.5;
nFFT = 2^ceil(log2(Nyquist));
%Delay = nFFT/(2*Fs)*10^3;
NBoot = 500; % number of voc ID permutation bootstraps for the significance of Info on coherence for each cell
% Lags = -Delay:Delay;
% Freqs = (0:ceil(length(Lags)/2)).* (2*Nyquist/length(Lags)); % Lags is a uneven number so F(i) = i*2*Nyquist/length(Lags)

fprintf(1, 'Running coherence4LMC on %s %s Self=%d\n', FeatureName, Session, Self)
CellTimer = tic();


%% Various checks (Number of vocalizations in the dataset.. etc)
if ~isfield(Cell, 'What')
    fprintf(1,'*** Problem with Cell: no what field!! ****\n')
    Coherence.Error = '*** No what field!! ****';
    fprintf(1, '------------------Done with Calculations %s %s Self=%d in %ds-------------------\n',  FeatureName, Session, Self,toc(CellTimer));
    return
else
    if Self % Select vocalizations from the requested session, emitted by self and that have the right delay before and after
        IndVoc = find(contains(Cell.What, 'Voc') .* contains(Cell.ExpType, Session(1)) .* contains(Cell.Who, 'self') .* (Cell.DelayBefore>=Delay) .* (Cell.DelayAfter>=Delay));
    else % Select vocalizations from the requested session, emitted by...
        % Others with no overlap with other vocalizations and that have
        % the right delay before and after, % ensure as well that there is
        % no noise on the microphone track (manual annotation) if we are
        % in free session, no check for operant
        IndVoc = find(contains(Cell.What, 'Voc') .* contains(Cell.ExpType, Session(1)) .*(~contains(Cell.Who, 'self')) .* (Cell.VocOverlap==0));
        if strcmp(Session(1), 'F')
            IndVoc = intersect(IndVoc, find(Cell.AudioQuality==1));
        end
    end
    if Trill % Only select Trill vocalizations
        IndVoc = intersect(IndVoc, find(contains(Cell.What, 'Tr')));
    end
    if Nvoc % We want to caclulate coherence on a fix number of vocalizations, randomly select these
        RandIndVoc = IndVoc(randperm(length(IndVoc)));
        IndVoc = RandIndVoc(1:Nvoc);
    end
    NStims = length(IndVoc);
    if NStims<10
        fprintf(1, '*** Problem with Cell: Not enough vocalizations only %d\n', NStims)
        Coherence.Error = sprintf('Not enough vocalizations only %d\n', NStims);
        fprintf(1, '------------------Done with Calculations %s %s Self=%d in %ds-------------------\n',  FeatureName, Session, Self,toc(CellTimer));
        return
    end
    StimDura = nan(NStims,1);
    for sss=1:NStims
        StimDura(sss) = round(length(Cell.BioSound{IndVoc(sss),2}.sound) ./(Cell.BioSound{IndVoc(sss),2}.samprate)*10^3);
    end
    if any(abs(StimDura - Cell.Duration(IndVoc))>1)
        fprintf(1,'*** Problem with Cell: duration inconcistency!! ****\n')
        Coherence.Error = 'Duration inconsistency';
        fprintf(1, '------------------Done with Calculations %s %s Self=%d in %ds-------------------\n',  FeatureName, Session, Self,toc(CellTimer));
        return
    end
end



%% Compute neural vectors
% Neural Data loop
% neural response is a vector that compile all spike counts starting
% at -200ms (Delay) before stim onset and stop at 200ms (Delay) after
% stim offset
[YPerStim, YPerStimt] = get_y_4Coherence(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),Delay,TR);
Y = [YPerStim{:}]';
if nansum(Y)<=2
    fprintf(1,'Cell: Non Spiking cell, No calculation!!\n')
    Coherence.Error = 'Non Spiking cell, No calculation!!';
    fprintf(1, '------------------Done with Calculations %s %s Self=%d in %ds-------------------\n',  FeatureName, Session, Self,toc(CellTimer));
    return
end

%% Calculate acoustic features input to the models
% acoustic data is a vector of the value of the acoustic feature sampled
% at 1/TRHz starting -200ms (Delay) before stim onset and stop at 200ms (Delay) after
% stim offset
DefaultVal = 0;%zero should be the default value for the amplitude, we know here that there is no sound
[XPerStim, XPerStimt] = get_x_4coherence(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Delay,TR,DefaultVal,FeatureName);
if Average
    Length = cellfun(@length, XPerStim);
    XPerStim_mat = zeros(NStims, max(Length));
    for Stim = 1:NStims
        XPerStim_mat(Stim,1:Length(Stim)) = XPerStim{Stim};
    end
    AvStim = nanmean(XPerStim_mat,1);
    for Stim = 1:NStims
        XPerStim{Stim} = AvStim(1:Length(Stim));
    end
end
X = [XPerStim{:}]';

if StimXYDataPlot
    if strcmp(FeatureName, 'amp') %#ok<UNRCH>
        F_high = 10000;
        ColBiosound = 2;
        FeatureNameLocal = 'Amplitude';
    elseif strcmp(FeatureName, 'SpectralMean')
        F_high = 50000;
        ColBiosound = 1;
        FeatureNameLocal = 'Spectral Mean';
    elseif strcmp(FeatureName, 'sal')
        F_high = 10000;
        ColBiosound = 2;
        FeatureNameLocal = 'saliency';
    end
    plotxyfeaturescoherence(Cell.BioSound(IndVoc,ColBiosound),YPerStim,YPerStimt, XPerStim,XPerStimt,TR,Delay,Cell.Duration(IndVoc),F_high,FeatureNameLocal)
end

%% Calculate coherence and coherency between the signals
[CoherencyT_unshifted, Freqs, CSR, CSR_up, CSR_low,~] = multitapercoherence_JN_fast([Y X],nFFT,Fs);

%% Calculate information on coherence and other parameters on coherence
Coherence=coherenceinfo_cal(CoherencyT_unshifted,Freqs, CSR, CSR_up, CSR_low, Nyquist,Fs,nFFT,TR,PlotCoherenceFig, FeatureName);
Coherence.Freqs =Freqs;
Coherence.LengthX = length(X);
Coherence.NStims = NStims;
if PlotCoherenceFig
    figure(3) %#ok<UNRCH>
%     suplabel(sprintf('Coherence Cell %d/%d', cc, NCells), 't');
    drawnow
end
%% Bootsrap the calculation of coherence and information with permutation...
% of vocalizations identity (X) BootstrapType==0
% or same vocalization ID but random onset time between 0 and BootstrapType
% ms or Bootstraptype(1) ms and BootstrapType(2)ms
if BootstrapType
    Info_boot = cell(1,NBoot);
    CohT_DelayAtzero_boot = cell(1,NBoot);
    CohT_WidthAtMaxPeak_boot = cell(1,NBoot);
    warning('off', 'signal:findpeaks:largeMinPeakHeight')
    if BootstrapType==4
        Yboot_All = get_y_4Coherence_rand(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),Delay,TR,length(Y), sum(Y),0, NBoot);
    end
    for bb=1:NBoot %parfor
        if ~mod(bb,NBoot/10)
            fprintf(1, 'Bootstrap %d  ', bb)
        end
        warning('off', 'signal:findpeaks:largeMinPeakHeight')
        if BootstrapType==1
            Xboot = [XPerStim{randperm(length(XPerStim))}]';
            Yboot = Y;
        elseif BootstrapType==2
            Yboot = YPerStim;
            for Stim=1:length(YPerStim)
                TimeShuffle = randi([1 length(YPerStim{Stim})]);
                Yboot{Stim} = [YPerStim{Stim}(TimeShuffle:end) YPerStim{Stim}(1:TimeShuffle-1)];  
            end
            Yboot = [Yboot{:}]';
            Xboot = X;
        elseif BootstrapType==3
            TimeShuffle = randi([1 length(Y)]);
            Yboot = [Y(TimeShuffle:end); Y(1:TimeShuffle-1)];  
            Xboot = X;
        elseif BootstrapType == 4
            Yboot = Yboot_All{bb}';
            Xboot = X;
        else
            error('Wrong input for BootstrapType')
        end
        [CoherencyT_unshifted, Freqs_b, CSR, ~, CSR_low,~] = multitapercoherence_JN_fast([Yboot Xboot],nFFT,Fs);
        [Bootstrap]=coherenceinfo_cal_bootstrap(CoherencyT_unshifted,Freqs_b, CSR, CSR_low, Nyquist,Fs,nFFT, TR);
        Info_boot{bb} = Bootstrap.Info;
        CohT_DelayAtzero_boot{bb} = Bootstrap.CoherencyT_DelayAtzero;
        CohT_WidthAtMaxPeak_boot{bb} = Bootstrap.CoherencyT_WidthAtMaxPeak;
    end
    warning('on', 'signal:findpeaks:largeMinPeakHeight')
    fprintf(1, 'DONE\n')
    if BootstrapType==1
        Coherence.Bootstrap.Info = [Info_boot{:}];
        Coherence.Bootstrap.CoherencyT_DelayAtzero = [CohT_DelayAtzero_boot{:}];
        Coherence.Bootstrap.CoherencyT_WidthAtMaxPeak = [CohT_WidthAtMaxPeak_boot{:}];
        Info_p = sum(([Info_boot{:}] - Coherence.Info)>= 0)/NBoot;
        Coherence.Info_p = Info_p;
    elseif BootstrapType==2
        Coherence.BootstrapTime.Info = [Info_boot{:}];
        Coherence.BootstrapTime.CoherencyT_DelayAtzero = [CohT_DelayAtzero_boot{:}];
        Coherence.BootstrapTime.CoherencyT_WidthAtMaxPeak = [CohT_WidthAtMaxPeak_boot{:}];
        Info_p = sum(([Info_boot{:}] - Coherence.Info)>= 0)/NBoot;
        Coherence.Info_pTime = Info_p;

    elseif BootstrapType==3
        Coherence.BootstrapFullTime.Info = [Info_boot{:}];
        Coherence.BootstrapFullTime.CoherencyT_DelayAtzero = [CohT_DelayAtzero_boot{:}];
        Coherence.BootstrapFullTime.CoherencyT_WidthAtMaxPeak = [CohT_WidthAtMaxPeak_boot{:}];
        Info_p = sum(([Info_boot{:}] - Coherence.Info)>= 0)/NBoot;
        Coherence.Info_pFullTime = Info_p;
        
    elseif BootstrapType==4
        Coherence.BootstrapRandSpikePerm.Info = [Info_boot{:}];
        Coherence.BootstrapRandSpikePerm.CoherencyT_DelayAtzero = [CohT_DelayAtzero_boot{:}];
        Coherence.BootstrapRandSpikePerm.CoherencyT_WidthAtMaxPeak = [CohT_WidthAtMaxPeak_boot{:}];
        Info_p = sum(([Info_boot{:}] - Coherence.Info)>= 0)/NBoot;
        Coherence.Info_pRandSpikePerm = Info_p;
    end

    if PlotCoherenceFig
        figure(4) %#ok<UNRCH>
        clf
        histogram([Info_boot{:}])
        hold on
        VL = vline(Coherence.Info, '-r');
        VL.LineWidth = 2;
        ylabel(sprintf('# bootstrap (total = %d)', NBoot))
        xlabel(sprintf('Information on coherence with sound %s', FeatureName))
        title(sprintf('Cell Significance of information with bootstraped permutation test p=%.2f', Info_p))
        drawnow
    end
end
%     keyboard
%     if Info(cc)>2
%         keyboard
%     end
fprintf(1, '------------------Done with Calculations %s %s Self=%d in %ds-------------------\n',  FeatureName, Session, Self,toc(CellTimer));



%% INTERNAL FUNCTIONS

    function [XPerStim, XPerStimt] = get_x_4coherence(BioSound, Duration, Delay,TR,DefaultVal,Feature, Overlap)
        DebugFig=0;
        if nargin<7
            Overlap = 0;
        end
        if length(Delay)==1
            Delay = [Delay Delay];
        end
        % calculate the cell array of acoustic data for each cell. For each
        % voc, XPerStim is a vector where each column correspond to the values of
        % the acoustic feature in the window Win starting Delay ms before time onset
        % of the vocalization.
        %  Fs = round(1/((TR-Overlap).*10^-3));
        % Vectors of the acoustic features
        XPerStim = cell(1,length(Duration));
        XPerStimt = cell(1,length(Duration));
        if ischar(DefaultVal) && strcmp(DefaultVal, 'mean')
            DefaultValue = 1;
        else
            DefaultValue = DefaultVal;
        end
        
        for stim = 1:length(Duration)
            if abs(round(length(BioSound{stim}.sound)/BioSound{stim}.samprate*10^3) - round(Duration(stim)))>1
                keyboard
            end
            % Get ready an output vector for the stim acoustic features that was sampled at 1000Hz
            XPerStim_temp = DefaultValue.*ones(1,round(1000 .* (Delay(1) + Duration(stim) + Delay(2)).*10^-3));
            FeatureVal = BioSound{stim}.(sprintf('%s',Feature));
            XPerStim_temp(Delay(1)+(1:length(FeatureVal)))=FeatureVal;
            if ~ischar(DefaultVal)
                XPerStim_temp(isnan(XPerStim_temp)) = DefaultValue;
            else
                XPerStim_temp(isnan(XPerStim_temp)) = nanmean(XPerStim_temp);
            end
            % Warning, resample gives very weird results here!!!
            %     [XPerStim{stim}] = resample(XPerStim_temp, Fs, 1000);
            %     [XPerStim{stim},XPerStimt{stim}] = resample(XPerStim_temp, (1:length(XPerStim_temp)),Fs*10^-3, 'spline');
            
            % I'm doing my own resampling
            % Time slots
            TimeBinsXOnsetInd = 1 :(TR-Overlap): (Delay(2) + Delay(1) + Duration(stim)); % These are slightly different than in get_Y_4GLM, because the times slot are used as indices in the vector and not as actuel time values!
            TimeBinsXOffsetInd = TimeBinsXOnsetInd + TR -1;
            TimeBinsXOnsetInd = TimeBinsXOnsetInd(TimeBinsXOffsetInd<=(Delay(2) + Delay(1) + Duration(stim))); % Only keep windows that are within the call
            TimeBinsXOffsetInd = TimeBinsXOffsetInd(TimeBinsXOffsetInd<=(Delay(2) + Delay(1) + Duration(stim))); % Only keep windows that are within the call
           
            TimeBinsXOnset = -Delay(1) :(TR-Overlap): (Delay(2) + Duration(stim));
            TimeBinsXOffset = TimeBinsXOnset + TR;
            TimeBinsXOnset = TimeBinsXOnset(TimeBinsXOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
            TimeBinsXOffset = TimeBinsXOffset(TimeBinsXOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
            XPerStimt{stim} = TimeBinsXOnset + (TimeBinsXOffset - TimeBinsXOnset)/2;
            
            XPerStim_resamp = nan(TR, length(TimeBinsXOnsetInd));
            for tt=1:length(TimeBinsXOnsetInd)
                XPerStim_resamp(:,tt) = XPerStim_temp(TimeBinsXOnsetInd(tt): TimeBinsXOffsetInd(tt))';
            end
            XPerStim{stim} = mean(XPerStim_resamp);
            
            
            
            if DebugFig
                figure(200) %#ok<UNRCH>
                clf
                plot((-Delay(1)+0.5):(Duration(stim)+Delay(2)),XPerStim_temp, 'LineWidth',2)
                xlabel('Time ms')
                ylabel(sprintf('%s',Feature))
                title(sprintf('Stim %d/%d',stim,length(Duration)));
                %         RemainTime = length(XPerStim_temp) - (length(XPerStim{stim})-1)*(1/Fs*10^3);
                %         XPerStimt{stim} = RemainTime/2+(1/Fs*10^3)*(0:(length(XPerStim{stim})-1));
                hold on
                plot(XPerStimt{stim},XPerStim{stim}, 'LineWidth',2)
                legend({'original' 'resampled'})
                pause(1)
            end
        end
        
    end



    function [YPerStim, YPerStimt] = get_y_4Coherence(SAT, Duration,Delay,TR,Overlap)
        DebugFig = 0;
        if nargin<5
            Overlap = 0;
        end
        if length(Delay)==1
            Delay = [Delay Delay];
        end
        % Calculate the time varying rate applying a gaussian window TR on the
        % spike pattern. The spike pattern considered starts -Delay ms
        % before the onset of the vocalization and stops Delay ms after the
        % offset of the vocalization
        YPerStim = cell(1,length(Duration));
        YPerStimt = cell(1,length(Duration));
        % Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
        nStd =(max(Duration) + Delay(1) + Delay(2))/10; % before set as 4
        Tau = (TR/2);
        T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
        Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
        Expwav = Expwav./sum(Expwav);
        % Frequency at which the neural data should be sampled
        FS = round(1/((TR-Overlap).*10^-3));
        % Loop through the stimuli and fill in the matrix
        for stim=1:length(Duration)
            % Time slots for the neural response
            TimeBinsY = -Delay(1) : (Delay(2) + Duration(stim));
            SpikePattern = zeros(1,length(TimeBinsY)-1);
            for isp = 1:length(SAT{stim})
                SpikeInd = round(SAT{stim}(isp));
                if (SpikeInd>=-Delay(1)) && (SpikeInd<(Delay(2) + Duration(stim)))
                    SpikePattern(SpikeInd + Delay(1) +1) = SpikePattern(SpikeInd + Delay(1) +1) +1;
                end
            end
            
            % Convolve with Gaussian to obtain our smooth time varying spike train
            % and resample if necessary
            if FS == 1000
                YPerStim{stim} = conv(SpikePattern, Expwav,'same');
                YPerStimt{stim} = TimeBinsY;
            else
                YPerStim_local = conv(SpikePattern, Expwav,'same');
                if sum(YPerStim_local)>0
                    YPerStim_local = YPerStim_local/sum(YPerStim_local)*sum(SpikePattern); % Make sure we keep the right number of sipkes after convolution!
                end
                
                % resampling function is really doing weird things at edges...
                % doing my own resampling
                TimeBinsYOnsetInd = 1 :(TR-Overlap): (Delay(2) + Delay(1) + Duration(stim)); % These are slightly different than in get_Y_4GLM, because the times slot are used as indices in the vector and not as actuel time values!
                TimeBinsYOffsetInd = TimeBinsYOnsetInd + TR -1;
                TimeBinsYOnsetInd = TimeBinsYOnsetInd(TimeBinsYOffsetInd<=(Delay(2) + Delay(1) + Duration(stim))); % Only keep windows that are within the call
                TimeBinsYOffsetInd = TimeBinsYOffsetInd(TimeBinsYOffsetInd<=(Delay(2) + Delay(1) + Duration(stim))); % Only keep windows that are within the call
                
                
                
                TimeBinsYOnset = -Delay(1) :(TR-Overlap): (Delay(2) + Duration(stim));
                TimeBinsYOffset = TimeBinsYOnset + TR;
                TimeBinsYOnset = TimeBinsYOnset(TimeBinsYOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
                TimeBinsYOffset = TimeBinsYOffset(TimeBinsYOffset<=(Delay(2) + Duration(stim))); % Only keep windows that are within the call
                YPerStimt{stim} = TimeBinsYOnset + (TimeBinsYOffset - TimeBinsYOnset)/2;
                
                YPerStim_resamp = nan(TR, length(TimeBinsYOnsetInd));
                for tt=1:length(TimeBinsYOnsetInd)
                    YPerStim_resamp(:,tt) = YPerStim_local(TimeBinsYOnsetInd(tt): TimeBinsYOffsetInd(tt))';
                end
                YPerStim{stim} = mean(YPerStim_resamp);
               
                
                if DebugFig
                    figure(200) %#ok<UNRCH>
                    clf
                    plot((-Delay(1)+0.5):(Duration(stim)+Delay(2)),YPerStim_local, 'LineWidth',2)
                    xlabel('Time ms')
                    ylabel('Spike Rate mHz (/ms)')
                    title(sprintf('Stim %d/%d',stim,length(Duration)));
                    %         RemainTime = length(XPerStim_temp) - (length(XPerStim{stim})-1)*(1/Fs*10^3);
                    %         XPerStimt{stim} = RemainTime/2+(1/Fs*10^3)*(0:(length(XPerStim{stim})-1));
                    hold on
                    plot(YPerStimt{stim},YPerStim{stim}, 'LineWidth',2)
                    legend({'original' 'resampled'}, 'AutoUpdate','off')
                    hold on
                    SpikeTimes = TimeBinsY(logical(SpikePattern));
                    for ss = 1:length(SpikeTimes)
                        V=vline(SpikeTimes(ss), 'k-');
                        V.LineWidth = 2;
                        hold on
                    end
                    pause(1)
                end
                
            end
            
            %     % change zero values for the smallest value under matlab.
            %     if sum(YPerStim{stim}==0)
            %         MinData = min(YPerStim{stim}(YPerStim{stim} ~=0));
            %         if ~isempty(MinData)
            %             YPerStim{stim}(YPerStim{stim}==0)=min(MinData,realmin('double'));
            %         else
            %             YPerStim{stim}(YPerStim{stim}==0)=realmin('double');
            %         end
            %     end
            
            % Make sure that the output mean(Y) = input mean(Y)
            %     if abs(round(sum(YPerStim{stim})*TR) - sum(SpikePattern))>TR/5
            %         warning('discrepancy in spike rate calculations larger than TR/5= %d?', TR/5)
            %         keyboard
            %     end
            
            %     if any(YPerStim{stim}<0)
            %         keyboard
            %     end
            
            
            %         % Time slots for the neural response
            %         TimeBinsY = -(Delay) : TR: (Delay + Duration(stim));
            %         YPerStim{stim} = nan(1,length(TimeBinsY)-1);
            %         for tt=1:(length(TimeBinsY)-1)
            %             % Find the number of spikes
            %             YPerStim{stim}(tt) = sum( (SAT{stim}>=TimeBinsY(tt)) .* (SAT{stim}<TimeBinsY(tt+1)));
            %         end
        end
        
    end



    function [Y_rand] = get_y_4Coherence_rand(SAT, Duration,Delay,TR,LengthY, NSpike,Overlap, NBoot)
        if nargin<7
            Overlap = 0;
        end
        if nargin<8
            NBoot = 500;
        end
        if length(Delay)==1
            Delay = [Delay Delay];
        end
        % Calculate the time varying rate for all concatenate stims, after
        % shuffling the spike accross the entire spike train, and applying
        % a gaussian window TR on the
        % spike pattern. The spike pattern considered starts -Delay ms
        % before the onset of the vocalization and stops Delay ms after the
        % offset of the vocalization
        % Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
        nStd =(max(Duration) + Delay(1) + Delay(2))/10; % before set as 4
        Tau = (TR/2);
        T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
        Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
        Expwav = Expwav./sum(Expwav);
        % Frequency at which the neural data should be sampled
        FS = round(1/((TR-Overlap).*10^-3));
        % Loop through the stimuli and fill in the matrix
        SpikePattern=cell(1, length(Duration));
        for stim=1:length(Duration)
            % Time slots for the neural response
            TimeBinsY = -Delay(1) : (Delay(2) + Duration(stim));
            SpikePattern{stim} = zeros(1,length(TimeBinsY)-1);
            for isp = 1:length(SAT{stim})
                SpikeInd = round(SAT{stim}(isp));
                if (SpikeInd>=-Delay(1)) && (SpikeInd<(Delay(2) + Duration(stim)))
                    SpikePattern{stim}(SpikeInd + Delay(1) +1) = SpikePattern{stim}(SpikeInd + Delay(1) +1) +1;
                end
            end
        end
        % Concatenate the spike trains in a single vector and suffle the
        % spike order
        SpikePattern = [SpikePattern{:}];
        
        Y_rand = cell(NBoot,1);
        parfor SP=1:NBoot
            SpikePattern_local = SpikePattern(randperm(length(SpikePattern)));
        
            % Convolve with Gaussian to obtain our smooth time varying spike train
            % and resample if necessary
            if FS == 1000
                Y_rand{SP} = conv(SpikePattern_local, Expwav,'same');
            else
                Y_rand_local = conv(SpikePattern_local, Expwav,'same');
                if sum(Y_rand_local)>0
                    Y_rand_local = Y_rand_local/sum(Y_rand_local)*sum(SpikePattern_local); % Make sure we keep the right number of sipkes after convolution!
                end

                % resampling function is really doing weird things at edges...
                % doing my own resampling
                TotalDuration = ((Delay(2) + Delay(1))*length(Duration) + sum(Duration));
                TimeBinsYOnsetInd = 1 :(TR-Overlap): TotalDuration; % These are slightly different than in get_Y_4GLM, because the times slot are used as indices in the vector and not as actuel time values!
                TimeBinsYOffsetInd = TimeBinsYOnsetInd + TR -1;
                TimeBinsYOnsetInd = TimeBinsYOnsetInd(TimeBinsYOffsetInd<=TotalDuration); % Only keep windows that are within the concatenated calls
                TimeBinsYOffsetInd = TimeBinsYOffsetInd(TimeBinsYOffsetInd<=TotalDuration); % Only keep windows that are within the concatenated calls


                Y_resamp = nan(TR, length(TimeBinsYOnsetInd));
                for tt=1:length(TimeBinsYOnsetInd)
                    Y_resamp(:,tt) = Y_rand_local(TimeBinsYOnsetInd(tt): TimeBinsYOffsetInd(tt))';
                end
                Y_resamp = mean(Y_resamp);
                Y_rand{SP} = Y_resamp(1:LengthY)./sum(Y_resamp(1:LengthY)).* NSpike;
            end
        end
        
        
    end



    function plotxyfeaturescoherence(BioSound,YPerStim,YPerStimt, XPerStim,XPerStimt,TR,Delay,Duration,F_high,FeatureName) %#ok<DEFNU>
        % This function plots for each stimulus the spectrogram with the
        % corresponding acoustic feature and spike rate
        if nargin<6
            FeatureName = 'Acoustic Feature';
        end
        if length(Delay)==1
            Delay = [Delay Delay];
        end
        DBNOISE =60;
        f_low = 0;
        for stim =1:length(BioSound)
            figure(1)
            clf
            ColorCode = get(groot,'DefaultAxesColorOrder');
            subplot(2,1,1)
            title(sprintf('vocalization %d/%d', stim,length(BioSound)))
            yyaxis left
            logB = BioSound{stim}.spectro;
            maxB = max(max(logB));
            minB = maxB-DBNOISE;
            imagesc(double(BioSound{stim}.to)*1000,double(BioSound{stim}.fo),logB);          % to is in seconds
            axis xy;
            caxis('manual');
            caxis([minB maxB]);
            cmap = spec_cmap();
            colormap(cmap);
            %         colorbar()
            v_axis = axis;
            v_axis(3)=f_low;
            v_axis(4)=F_high;
            axis(v_axis);
            xlabel('time (ms)'), ylabel('Frequency');
            
            
            hold on
            yyaxis right
            plot(XPerStimt{stim}, XPerStim{stim}, 'Color',ColorCode(2,:),'LineStyle','-', 'Marker','*', 'LineWidth',2)
            ylabel(sprintf('%s', FeatureName))
            if strcmp(FeatureName,'Spectral Mean')
                ylim(v_axis(3:4))
            elseif strcmp(FeatureName,'Saliency')
                YLim = get(gca, 'YLim');
                ylim([0 YLim(2)])
            elseif strcmp(FeatureName,'Amplitude')
                YLim = get(gca, 'YLim');
                ylim([0 YLim(2)])
            end
            
            xlim([-Delay(1) Duration(stim) + Delay(2)])
            hold off
            
            subplot(2,1,2)
            yyaxis left
            if length(YPerStim{stim})<4
                bar(YPerStimt{stim},YPerStim{stim}*10^3, 'BarWidth',0.2)
            else
                bar(YPerStimt{stim},YPerStim{stim}*10^3)
            end
            ylabel(sprintf('Spikes/s or Hz (%d ms Gauss win)', TR))
            YLIM = [0 500];
            ylim(YLIM)
            xlabel('Time (ms)')
            hold on
            line([0 Duration(stim)], [0.975 0.975].*YLIM(2), 'Color','k', 'LineWidth',12)
            hold on
            text(Duration(stim)/2.3,0.975.*YLIM(2),'Vocalization', 'Color', [1 1 1])
            hold on
            yyaxis right
            plot(XPerStimt{stim}, XPerStim{stim}, 'Color',ColorCode(2,:),'LineStyle','-', 'Marker','*', 'LineWidth',2)
            %     plot(TR/2 + (0:TR:((IndMax-1)*TR)), StimFeature(1:IndMax), 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
            ylabel(sprintf('%s', FeatureName))
            if strcmp(FeatureName,'Spectral Mean')
                ylim(v_axis(3:4))
            elseif strcmp(FeatureName,'Saliency')
                YLim = get(gca, 'YLim');
                ylim([0 YLim(2)])
            elseif strcmp(FeatureName,'Amplitude')
                YLim = get(gca, 'YLim');
                ylim([0 YLim(2)])
            end
            xlim([-Delay(1) Duration(stim) + Delay(2)])
            hold off
            pause(1)
        end
    end


end
