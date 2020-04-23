addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'));
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'));
%% Data Info
Filename = '190129_1023_VocExtractData_200.mat';
Path = '/Users/elie/Documents/ManipBats/LMC/190110_59882_11689_HoHa/20190129_data';
D1=load(fullfile(Path,[Filename(1:(end-8)) '.mat']));
D2 = load(fullfile(Path,Filename));

% Cell #
CellNum = 5;

% Delay inherited from the extraction and the merge threshold of
% vocalization in who_calls.m
Delay = str2double(Filename(end-6 : end-4));

% Number of stimuli in the dataset
NStims = length(D2.SpikeTimesVoc.Logger16.VocDuration);
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
Win = 5; % time in ms

%% Acoustic feature parameters
% Bandpass filter for the enveloppe calculation
BandPassFilter = [200 90000];
Fhigh_powers =[20 50 100]; % Frequency upper bound for calculating the envelope (time running RMS)
Fs_env = 1/Win*10^3; % Sample frequency of the enveloppe such as to have one value per Win
% design bandpass filter of raw ambient recording
[z,p,k] = butter(6,BandPassFilter/(D1.FS/2),'bandpass');
sos_raw_band = zp2sos(z,p,k);

%% Variable initialization
% organize acoustic data as a matrix where each
% column corresponds to amp env at t-100:Win:t+100, t being
% time of neural window;
% neural response is a vector that compile all spike counts starting
% at -100ms before stim onset and stop at 100ms after
% stim offset
% Calculate the number of rows of the matrix and  vector for preallocation of space
TotalDataTime = sum(D2.SpikeTimesVoc.Logger16.VocDuration) + Delay*NStims;
Nrows = round(TotalDataTime/Win);

% Time slots of the acoustic features
TimeBinsX = -(Delay/2) : Win : (Delay/2);
NTimeBinsX =length(TimeBinsX);
X_Amp = cell(length(Fhigh_powers),1);
X_Pow = cell(length(Fhigh_powers),1);
Y = nan(Nrows,1);
Amp_env_Mic = cell(NStims,1);
Power_env_Mic = cell(NStims,1);
YLocal = cell(NStims,1);
for ff=1:length(Fhigh_powers)
    X_Amp{ff} = nan(Nrows, NTimeBinsX);
    X_Pow{ff} = nan(Nrows, NTimeBinsX);
end

%% Gather the data
ii = 0; % This count the number of rows in the matrix and vector (number of samples)
% Loop through the stimuli and fill in the matrix
for stim=1:NStims
    % Calculate the acoustic feature
    % Temporal Enveloppe
    % bandpass filter the ambient mic recording
    Filt_RawVoc = filtfilt(sos_raw_band,1,D2.SpikeTimesVoc.Logger16.VocWave{stim});
    Amp_env_Mic{stim} = cell(length(Fhigh_powers),1);
    Power_env_Mic{stim} = cell(length(Fhigh_powers),1);
    for ff=1:length(Fhigh_powers)
        [Amp_env_Mic{stim}{ff}, Power_env_Mic{stim}{ff}] = running_rms(Filt_RawVoc, D1.FS, Fhigh_powers(ff), Fs_env);
    end
    
    % Time slots for the neural response
    TimeBinsY = -(Delay/2) : Win : (Delay/2 + D2.SpikeTimesVoc.Logger16.VocDuration(stim));
    YLocal{stim} = nan(length(TimeBinsY)-1,1);
    
    for tt=1:(length(TimeBinsY)-1)
        ii = ii+1;
        % Find the number of spikes
        YLocal{stim}(tt) = sum( (TimeBinsY(tt)<=D2.SpikeTimesVoc.Logger16.SpikesTimes_VocCall{stim,CellNum}) .* (TimeBinsY(tt+1)>D2.SpikeTimesVoc.Logger16.SpikesTimes_VocCall{stim,CellNum}));
        Y(ii) = YLocal{stim}(tt);
        % enter the corresponding values for the envelope
        for ff=1:length(Fhigh_powers)
            X_Amp{ff}(ii,:) = Amp_env_Mic{stim}{ff}((1:NTimeBinsX) + tt-1);
            X_Pow{ff}(ii,:) = Power_env_Mic{stim}{ff}((1:NTimeBinsX) + tt-1);
        end
    end
    close all
    figure(1)
    ColorCode = get(groot,'DefaultAxesColorOrder');
    clf
    subplot(2,1,1)
    yyaxis right
    for ff=1:length(Fhigh_powers)
        plot(Amp_env_Mic{stim}{ff}((Delay/2/Win):(end-Delay/2/Win)), 'Color',ColorCode(ff+1,:),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('FHigh %d Hz',Fhigh_powers(ff)))
        hold on
    end
    legend('show')
    legend('AutoUpdate', 'off')
    ylabel('Sound RMS')
    yyaxis left
    bar(YLocal{stim})
    ylabel(sprintf('Number of spikes per %d ms bin', Win))
    ylim([0 Win])
    XTick = get(gca, 'XTickLabel');
    XTick = cellfun(@str2double, XTick) * Win;
    set(gca,'XTickLabel',XTick)
    xlabel('Time (ms)')
    title(sprintf('Stim %d/%d',stim,NStims))
    
    subplot(2,1,2)
    yyaxis right
    for ff=1:length(Fhigh_powers)
        plot(Power_env_Mic{stim}{ff}((Delay/2/Win):(end-Delay/2/Win)), 'Color',ColorCode(ff+1,:),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('FHigh %d Hz',Fhigh_powers(ff)))
        hold on
    end
    legend('show')
    legend('AutoUpdate', 'off')
    ylabel('Sound Power')
    yyaxis left
    bar(YLocal{stim})
    ylabel(sprintf('Number of spikes per %d ms bin', Win))
    ylim([0 Win])
    XTick = get(gca, 'XTickLabel');
    XTick = cellfun(@str2double, XTick) * Win;
    set(gca,'XTickLabel',XTick)
    xlabel('Time (ms)')
    hold off
    pause(1)
end

%% Calculate the model
% get rid of Nans
Nan_Ind=find(isnan(Y));
Y(Nan_Ind) = [];
for ff=1:length(Fhigh_powers)
    X_Amp{ff}(Nan_Ind,:)=[];
    X_Pow{ff}(Nan_Ind,:)=[];
end

%% Run ridge GLM Poisson on power
B_Env = cell(length(Fhigh_powers),1);
Dev_Env = nan(length(Fhigh_powers),1);
FitInfo_Env = cell(length(Fhigh_powers),1);
figure()
for ff=1:length(Fhigh_powers)
%     [B_Env{ff}, Dev_Env(ff)]=glmfit(X_Pow{ff},Y,ParamModel.DISTR,'link',ParamModel.LINK);
    [B_Env{ff}, FitInfo_Env{ff}]=lassoglm(X_Pow{ff},Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
    % find the model with the minimum of deviance (best lambda)
    [BestDev,BestMod] = min(FitInfo_Env{ff}.Deviance);
    plot(TimeBinsX,B_Env{ff}(:,BestMod),'Color',ColorCode(ff+1,:),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('FHigh %d Hz Dev = %.1f',Fhigh_powers(ff),BestDev))
    hold on
end
legend('show')
% XTick = get(gca, 'XTickLabel');
% XTick = cellfun(@str2double, XTick) * Win;
% set(gca,'XTickLabel',XTick)
xlabel('Time (ms)')
title('Poisson Ridge regression on Power')

%% Run ridge GLM Poisson on change of power
B_Env = cell(length(Fhigh_powers),1);
Dev_Env = nan(length(Fhigh_powers),1);
FitInfo_Env = cell(length(Fhigh_powers),1);
figure()
for ff=1:length(Fhigh_powers)
    X_local = abs(diff(X_Pow{ff}')');
%     [B_Env{ff}, Dev_Env(ff)]=glmfit(X_Pow{ff},Y,ParamModel.DISTR,'link',ParamModel.LINK);
    [B_Env{ff}, FitInfo_Env{ff}]=lassoglm(X_local,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
    % find the model with the minimum of deviance (best lambda)
    [BestDev,BestMod] = min(FitInfo_Env{ff}.Deviance);
    plot(TimeBinsX(2:end),B_Env{ff}(:,BestMod),'Color',ColorCode(ff+1,:),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('FHigh %d Hz Dev = %.1f',Fhigh_powers(ff),BestDev))
    hold on
end
legend('show')
% XTick = get(gca, 'XTickLabel');
% XTick = cellfun(@str2double, XTick) * Win;
% set(gca,'XTickLabel',XTick)
xlabel('Time (ms)')
title('Poisson Ridge regression on change of Power')


%% run logistic model
Y01 = Y>0;
B_Env01 = cell(length(Fhigh_powers),1);
Dev_Env01 = nan(length(Fhigh_powers),1);
FitInfo_Env01 = cell(length(Fhigh_powers),1);
figure()
for ff=1:length(Fhigh_powers)
%     [B_Env01{ff}, Dev_Env01(ff)]=glmfit(X_Pow{ff},Y01,'binomial','link','logit');
    [B_Env01{ff}, FitInfo_Env01{ff}]=lassoglm(X_Pow{ff},Y01,'binomial','Alpha', ParamModel.Alpha,'Link','logit','NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
    % find the model with the minimum of deviance (best lambda)
    [BestDev,BestMod] = min(FitInfo_Env01{ff}.Deviance);
    plot(TimeBinsX,B_Env{ff}(:,BestMod),'Color',ColorCode(ff+1,:),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('FHigh %d Hz Dev = %.1f',Fhigh_powers(ff),BestDev))
    hold on
end
legend('show')
XTick = get(gca, 'XTickLabel');
XTick = cellfun(@str2double, XTick) * Win;
set(gca,'XTickLabel',XTick)
xlabel('Time (ms)')
title('Logistic Ridge regression on Power')
%     [B_Env{ff}, FitInfo_Env]=lassoglm(X,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);

% lassoPlot(B_Env,FitInfo_Env,'PlotType','Lambda','XScale','log')
    
    
    

