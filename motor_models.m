addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'));
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'));
DatFig=0; %Set to 1 to see input data figures
%% Data Info
Filename = '59834_20190611_SSS_1-97.mat';
% Path = '/Users/elie/Documents/LMCResults/';
Path = '/Users/elie/Documents/ManipBats/LMC/ResultsFiles/';
load(fullfile(Path,Filename));


% Number of vocalizations in the dataset
IndVoc = find(contains(What, 'Voc') .* (Duration>20)); % No saliency calculated when duration <20ms
NStims = length(IndVoc);
%% Model parameters
% Assumption of stationarity over time
ParamModel.LAMBDARATIO=1e-4;
ParamModel.NUMLAMBDA=10;%25?
ParamModel.LINK='log';
ParamModel.DISTR='poisson';
% Determine a list of alpha (parameter that range the regularization
% betweeen ridge (L2, alpha=0) and lasso (L1, alpha =1))
ParamModel.Alpha=0.001; % STRFs are easier to interpret using ridge than using Lasso and deviances are similar.

% Define the time resolution at which analysis should be done
TRs = [1 5 10]; % time in ms

Win =200;%value in ms for the duration of the snippet of sound which
% acoustic features are used to predict the neural response at a given time
% t. This will be the size of the xaxis of the MRF.

Delays = [0 5 10 25 50 75 100];% ms
%The neural activity at time t is predicted by a time varying
% acoustic feature that starts at t+Delay ms

% %% Acoustic feature parameters
% % Bandpass filter for the enveloppe calculation
% BandPassFilter = [200 90000];
% Fhigh_powers =[20 50 100]; % Frequency upper bound for calculating the envelope (time running RMS)
% Fs_env = 1/Win*10^3; % Sample frequency of the enveloppe such as to have one value per Win
% % design bandpass filter of raw ambient recording
% [z,p,k] = butter(6,BandPassFilter/(D1.FS/2),'bandpass');
% sos_raw_band = zp2sos(z,p,k);

%% Variable initialization
% organize acoustic data as a matrix where each
% column corresponds to amp env at t-100:Win:t+100, t being
% time of neural window;
% neural response is a vector that compile all spike counts starting
% at -100ms before stim onset and stop at 100ms after
% stim offset
BestDevSpecMean = nan(length(TRs),length(Delays)); % contains the deviance of the spectral model
BestDevSal = nan(length(TRs),length(Delays));  % contains the deviance of the saliency model
% BestDevSpecMed = nan(length(TRs),length(Delays)); % contains the deviance of the spectral model

for tr = 1:length(TRs)
    TR = TRs(tr);
    fprintf(1,'Models with Time resolution %d ms (%d/%d)\n', TR, tr, length(TRs));
    for dd = 1:length(Delays)
        Delay = Delays(dd);
        fprintf(1,'Models with Delay= %d ms (%d/%d)\n', Delay,dd,length(Delays));
        %% Gather the data

        % Neural Data loop
        YPerStim = get_y(SpikesArrivalTimes_Behav(IndVoc), Duration(IndVoc),Delay,TR);
        Y = [YPerStim{:}]'; 

        % Acoustic feature loop, first column of biosound is microphone second is
        % piezo
        XSpecMeanPerStim = get_x(BioSound(IndVoc,1), Duration(IndVoc), Win, TR, Delay,'SpectralMean');
        XSpecMean = [XSpecMeanPerStim{:}]';
        XSaliencyPerStim = get_x(BioSound(IndVoc,2), Duration(IndVoc), Win, TR, Delay,'sal');
        XSaliency = [XSaliencyPerStim{:}]';
%         XSpecMedPerStim = get_x(BioSound(IndVoc,1), Duration(IndVoc), Win, TR, Delay,'Q2t');
%         XSpecMed = [XSpecMedPerStim{:}]';

        %% Plot of features and neuronal response if requested (BioSound,YPerStim,XPerStim, TR,Delay,F_high,FeatureName)
        if DatFig
            plotxyfeatures(BioSound(IndVoc,1),YPerStim,XSpecMeanPerStim,TR,Delay,Duration(IndVoc), 50000,'Spectral Mean')
            plotxyfeatures(BioSound(IndVoc,2),YPerStim,XSaliencyPerStim,TR,Delay,Duration(IndVoc), 10000,'Saliency')
%         plotxyfeatures(BioSound(IndVoc,1),YPerStim,XSpecMedPerStim,TR,Delay,Duration(IndVoc), 50000,'Spectral Median')
        end
        
        %% Calculate the model
        % get rid of Nans
        Nan_Ind=find(isnan(Y));
        Y(Nan_Ind) = [];
        XSpecMean(Nan_Ind,:) =[];
        XSaliency(Nan_Ind,:) =[];
%         XSpecMed(Nan_Ind,:) =[];
        
        %% Run ridge GLM Poisson on acoustic features

        % spectral mean predicting Y
        [BSpecMean, FitInfo_SpecMean]=lassoglm(XSpecMean,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevSpecMean(tr,dd),BestModSpecMean] = min(FitInfo_SpecMean.Deviance);

        % pitch saliency predicting Y
        [BSal, FitInfo_Sal]=lassoglm(XSaliency,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevSal(tr,dd),BestModSal] = min(FitInfo_Sal.Deviance);
        
%         % spectral median predicting Y
%         [BSpecMed, FitInfo_SpecMed]=lassoglm(XSpecMed,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
%         % find the model with the minimum of deviance (best lambda)
%         [BestDevSpecMed(tr,dd),BestModSpecMed] = min(FitInfo_SpecMean.Deviance);

        TimeBinsX = 0 : TR : Win;
        figure()
        plot(TimeBinsX,BSpecMean(:,BestModSpecMean),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Spectral Mean Dev = %.1f',BestDevSpecMean(tr,dd)))
        hold on
        plot(TimeBinsX,BSal(:,BestModSal),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Saliency Dev = %.1f',BestDevSal(tr,dd)))
        hold on
%         plot(TimeBinsX,BSpecMed(:,BestModSpecMed),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Spectral Median Dev = %.1f',BestDevSpecMed(tr,dd)))

        legend('show')
        % XTick = get(gca, 'XTickLabel');
        % XTick = cellfun(@str2double, XTick) * Win;
        % set(gca,'XTickLabel',XTick)
        xlabel('Time (ms)')
        title(sprintf('Poisson Ridge regression on Acoustic Features Delay = %dms Win=%dms', Delay, TR))
        hold off
        pause(1)
    end
end

figure()
imagesc(Delays,TRs,BestDevSpecMean)
xlabel('Delays to voc onset in ms')
ylabel('Time resolution in ms')
title('Deviance spectral mean model')
colorbar()

figure()
imagesc(Delays,TRs,BestDevSal)
xlabel('Delays to voc onset in ms')
ylabel('Time resolution in ms')
title('Deviance pitch saliency model')
colorbar()

% figure()
% imagesc(Delays,TRs,BestDevSpecMed)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance spectral median model')

% %% Run ridge GLM Poisson on change of power
% B = cell(length(Fhigh_powers),1);
% Dev_Env = nan(length(Fhigh_powers),1);
% FitInfo_Env = cell(length(Fhigh_powers),1);
% figure()
% for ff=1:length(Fhigh_powers)
%     X_local = abs(diff(X_Pow{ff}')');
% %     [B_Env{ff}, Dev_Env(ff)]=glmfit(X_Pow{ff},Y,ParamModel.DISTR,'link',ParamModel.LINK);
%     [B{ff}, FitInfo_Env{ff}]=lassoglm(X_local,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
%     % find the model with the minimum of deviance (best lambda)
%     [BestDevSal,BestModSal] = min(FitInfo_Env{ff}.Deviance);
%     plot(TimeBinsX(2:end),B{ff}(:,BestModSal),'Color',ColorCode(ff+1,:),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('FHigh %d Hz Dev = %.1f',Fhigh_powers(ff),BestDevSal))
%     hold on
% end
% legend('show')
% % XTick = get(gca, 'XTickLabel');
% % XTick = cellfun(@str2double, XTick) * Win;
% % set(gca,'XTickLabel',XTick)
% xlabel('Time (ms)')
% title('Poisson Ridge regression on change of Power')
% 
% 
% %% run logistic model
% Y01 = Y>0;
% B_Env01 = cell(length(Fhigh_powers),1);
% Dev_Env01 = nan(length(Fhigh_powers),1);
% FitInfo_Env01 = cell(length(Fhigh_powers),1);
% figure()
% for ff=1:length(Fhigh_powers)
% %     [B_Env01{ff}, Dev_Env01(ff)]=glmfit(X_Pow{ff},Y01,'binomial','link','logit');
%     [B_Env01{ff}, FitInfo_Env01{ff}]=lassoglm(X_Pow{ff},Y01,'binomial','Alpha', ParamModel.Alpha,'Link','logit','NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
%     % find the model with the minimum of deviance (best lambda)
%     [BestDevSal,BestModSal] = min(FitInfo_Env01{ff}.Deviance);
%     plot(TimeBinsX,B{ff}(:,BestModSal),'Color',ColorCode(ff+1,:),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('FHigh %d Hz Dev = %.1f',Fhigh_powers(ff),BestDevSal))
%     hold on
% end
% legend('show')
% XTick = get(gca, 'XTickLabel');
% XTick = cellfun(@str2double, XTick) * TR;
% set(gca,'XTickLabel',XTick)
% xlabel('Time (ms)')
% title('Logistic Ridge regression on Power')
% %     [B_Env{ff}, FitInfo_Env]=lassoglm(X,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
% 
% % lassoPlot(B_Env,FitInfo_Env,'PlotType','Lambda','XScale','log')
%     
    
   

%% INTERNAL FUNCTIONS

function [XPerStim] = get_x(BioSound, Duration, Win, TR,Delay,Feature)
    % Time slots of the acoustic features
    TimeBinsX = 0 : TR : Win;
    NTimeBinsX =length(TimeBinsX);
    XPerStim = cell(1,length(Duration));
    for stim = 1:length(Duration)
        % Time slots of the corresponding neural response
        TimeBinsY = -(Delay) : TR : (Delay + Duration(stim));
        
        % Get ready the stim acoustic features that was sampled at 1000Hz
        FeatureVal = BioSound{stim}.(sprintf('%s',Feature));
        FeatureValResampled = resample(FeatureVal,1,TR); % Win in ms so resampling from 1 value per ms to 1/Win value per ms or 1 value per win ms
        ZPaddFeature = zeros(1,NTimeBinsX+length(TimeBinsY)-1);
        ZPaddFeature(1:length(FeatureValResampled))=FeatureValResampled;
        
        XPerStim{stim} = zeros(length(TimeBinsX),length(TimeBinsY)-1);
        for tt=1:(length(TimeBinsY)-1)
            XPerStim{stim}(:,tt) = ZPaddFeature((1:NTimeBinsX) + tt-1);
        end
    end

end

function [YPerStim] = get_y(SAT, Duration, Delay,TR)
    YPerStim = cell(1,length(Duration));
    % Loop through the stimuli and fill in the matrix
    for stim=1:length(Duration)
        % Time slots for the neural response
        TimeBinsY = -(Delay) : TR : (Delay + Duration(stim));
        YPerStim{stim} = nan(1,length(TimeBinsY)-1);
    
        for tt=1:(length(TimeBinsY)-1)
            % Find the number of spikes
            YPerStim{stim}(tt) = sum( (SAT{stim}>=TimeBinsY(tt)) .* (SAT{stim}<TimeBinsY(tt+1)));
        end
    end

end

%     % Calculate the acoustic feature
%     % Temporal Enveloppe
%     % bandpass filter the ambient mic recording
%     Filt_RawVoc = filtfilt(sos_raw_band,1,D2.SpikeTimesVoc.Logger16.VocWave{stim});
%     Amp_env_Mic{stim} = cell(length(Fhigh_powers),1);
%     Power_env_Mic{stim} = cell(length(Fhigh_powers),1);
%     for ff=1:length(Fhigh_powers)
%         [Amp_env_Mic{stim}{ff}, Power_env_Mic{stim}{ff}] = running_rms(Filt_RawVoc, D1.FS, Fhigh_powers(ff), Fs_env);
%     end
 
function plotxyfeatures(BioSound,YPerStim,XPerStim, TR,Delay,Duration,F_high,FeatureName)
if nargin<6
    FeatureName = 'Acoustic Feature';
end
DBNOISE =12;
f_low = 0;
close all
for stim =1:length(BioSound)
    figure(1)
    clf
    ColorCode = get(groot,'DefaultAxesColorOrder');
    subplot(2,1,1)
    title(sprintf('vocalization %d/%d', stim,length(BioSound)))
    yyaxis left
    logB = - 20*log10(abs(double(BioSound{stim}.spectro)));
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
    StimFeature = [XPerStim{stim}(1,:) XPerStim{stim}(2:end,end)'];
    IndMax = floor(max(double(BioSound{stim}.to)*1000)/TR);
    plot(TR/2 + (0:TR:((IndMax-1)*TR)), StimFeature(1:IndMax), 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    ylabel(sprintf('%s', FeatureName))
    if strcmp(FeatureName,'Spectral Mean')
        ylim(v_axis(3:4))
    end
    hold off
    
    subplot(2,1,2)
    yyaxis left
    bar(TR/2+(-(Delay) : TR : (Delay + Duration(stim)-TR)),YPerStim{stim})
    ylabel(sprintf('Number of spikes per %d ms bin', TR))
    ylim([0 TR])
    xlabel('Time (ms)')
    hold on
    yyaxis right
    plot(TR/2 + (0:TR:((IndMax-1)*TR)), StimFeature(1:IndMax), 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    ylabel(sprintf('%s', FeatureName))
    XLIM = get(gca,'XLim');
    hold off
    subplot(2,1,1)
    set(gca,'XLim',XLIM)
    pause(1)
end
end

 

