addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'));
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'));
DatFig=0; %Set to 1 to see input data figures

%% Listing datacells
%Filename = '59834_20190611_SSS_1-97.mat';
% Filename = '59834_20190610_SSS_1-130.mat';
Path = '/Users/elie/Documents/LMCResults/';
% Path = '/Users/elie/Documents/ManipBats/LMC/ResultsFiles/';

AllFiles = dir(fullfile(Path,'*.mat'));
Files2run = zeros(length(AllFiles),1);
for ff=1:length(AllFiles)
    if length(strfind(AllFiles(ff).name, '_'))==3
        Files2run(ff) = 1;
    end
end
CellsPath = AllFiles(logical(Files2run));
%% Model parameters
% Assumption of stationarity over time
ParamModel.LAMBDARATIO=1e-4;
ParamModel.NUMLAMBDA=10;%25?
ParamModel.LINK='log';
ParamModel.DISTR='poisson';
% Determine a list of alpha (parameter that range the regularization
% betweeen ridge (L2, alpha=0) and lasso (L1, alpha =1))
ParamModel.Alpha=0.001; % STRFs are easier to interpret using ridge than using Lasso and deviances are similar.

% Define the time resolution at which neural density estimates should be calculated
% TRs = [1 5 10 15 20 25 30 35 40 45 50]; % time in ms
% Based on first results of calculations for a subset of cells, I decide to
% rping down the number of time resolution to explore
TRs = [1 5 10 15 20]; % time in ms

Win =300;%value in ms for the duration of the snippet of sound which
% acoustic features are used to predict the neural response at a given time
% t. This will be the size of the xaxis of the MRF.

Delay = 100;% ms
%The neural activity at time t is predicted by a time varying
% acoustic feature that starts at t-Delay ms


%% Running through cells
NCells = length(CellsPath);
%MotorModels = cell(NCells,1);
load(fullfile(Path,'MotorModelsAllCells'),'MotorModels','CellsPath');
parfor cc=1:NCells
    if ~isempty(MotorModels{cc}) % This one was already calculated
        fprintf(1, 'Cell %d/%d Already calculated\n',cc,NCells)
        continue
    end
    Cell = load(fullfile(CellsPath(cc).folder,CellsPath(cc).name));
    
    
    % Number of vocalizations in the dataset
    if ~isfield(Cell, 'What')
        fprintf(1,'*** . Problem with Cell %d, no what field!! ****\n', cc)
    end
    IndVoc = find(contains(Cell.What, 'Voc') .* (Cell.Duration>20)); % No saliency calculated when duration <20ms
    NStims = length(IndVoc);
    
    %% Variable organization and running models
    % organize acoustic data as a matrix where each
    % column corresponds to amp env at t-100:Win:t+100, t being
    % time of neural window;
    % neural response is a vector that compile all spike counts starting
    % at -100ms before stim onset and stop at 100ms after
    % stim offset
    BestDevSpecMean = nan(length(TRs),1); % contains the deviance of the spectral mean model
    BestDevSal = nan(length(TRs),1);  % contains the deviance of the saliency model
    BestDevAmp = nan(length(TRs),1);  % contains the deviance of the amplitude model
    % BestDevSpecMed = nan(length(TRs),1); % contains the deviance of the spectral median model
    BestDevNull = nan(length(TRs),1);  % contains the deviance of the null model
    
    BSpecMean = nan(length(TRs),2*Win+1); % contains the Betas of the spectral mean model with the best deviance
    BSal = nan(length(TRs),2*Win+1);  % contains the Betas of the saliency model with the best deviance
    BAmp = nan(length(TRs),Win+1);  % contains the Betas of the amplitude model with the best deviance
    % BestDevSpecMed = nan(length(TRs),1); % contains the Betas of the spectral
    % median model with the best deviance
    % BNull = nan(length(TRs),Win+1);  % contains the Beta of the null model
    BNull = nan(length(TRs),1);  % contains the Beta of the null model
    TicToc = nan(length(TRs),1);  % duration of each loop of models
    MeanY = nan(length(TRs),1);
    
    for tr = 1:length(TRs)
        TR = TRs(tr);
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d)\n', cc,NCells,TR, tr, length(TRs));
        %% Gather the data
        TimerLoop =tic;
        % Neural Data loop
        YPerStim = get_y(Cell.SpikesArrivalTimes_Behav(IndVoc), Cell.Duration(IndVoc),Win,Delay,TR);
        Y = [YPerStim{:}]';
        
        MeanY(tr)=nanmean(Y);
        % check that the neural response is high enough
        if nanmean(Y)<10^-2
            fprintf(1, 'Cell %d/%d Models with Time resolution %d ms (%d/%d) -> Rate is too low, not calculating models to avoid convergence issues\n',cc,NCells,TR, tr, length(TRs))
            continue
        end
        
        
        % Acoustic feature loop, first column of biosound is microphone second is
        % piezo
        XSpecMeanPerStim = get_x(Cell.BioSound(IndVoc,1), Cell.Duration(IndVoc), Win, Delay,'SpectralMean');
        XSpecMean = [XSpecMeanPerStim{:}]';
        XSpecMean(isnan(XSpecMean))=nanmean(reshape(XSpecMean,numel(XSpecMean),1)); % Change the padding with NaN to the mean for the SpecMean for silence.
        XSaliencyPerStim = get_x(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Win, Delay,'sal');
        XSaliency = [XSaliencyPerStim{:}]';
        XSaliency(isnan(XSaliency))=0; % Change the padding with NaN to zeros for the Saliency, we know Silence is not harmonic
        XAmpPerStim = get_x(Cell.BioSound(IndVoc,2), Cell.Duration(IndVoc), Win, Delay,'amp');
        XAmp = [XAmpPerStim{:}]';
        XAmp(isnan(XAmp))=0; % Change the padding with NaN to zeros for the amplitude, we know here that there is no sound
        %         XSpecMedPerStim = get_x(BioSound(IndVoc,1), Duration(IndVoc), Win, TR, Delay,'Q2t');
        %         XSpecMed = [XSpecMedPerStim{:}]';
        
        %% Plot of features and neuronal response if requested (BioSound,YPerStim,XPerStim, TR,Delay,F_high,FeatureName)
        if DatFig
            plotxyfeatures(BioSound(IndVoc,1),YPerStim,XSpecMeanPerStim,TR,Win,Delay,Duration(IndVoc), 50000,'Spectral Mean')
            plotxyfeatures(BioSound(IndVoc,2),YPerStim,XSaliencyPerStim,TR,Win,Delay,Duration(IndVoc), 10000,'Saliency')
            plotxyfeatures(BioSound(IndVoc,2),YPerStim,XAmpPerStim,TR,Win,Delay,Duration(IndVoc), 10000,'Amp')
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
        [BSpecMean_local, FitInfo_SpecMean]=lassoglm([XAmp XSpecMean],Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',1,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevSpecMean(tr),BestModSpecMean] = min(FitInfo_SpecMean.Deviance);
        BSpecMean(tr,2:end) = BSpecMean_local(:,BestModSpecMean);
        BSpecMean(tr,1) = FitInfo_SpecMean.Intercept(BestModSpecMean);
        
        % pitch saliency predicting Y
        [BSal_local, FitInfo_Sal]=lassoglm([XAmp XSaliency],Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',1,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevSal(tr),BestModSal] = min(FitInfo_Sal.Deviance);
        BSal(tr,2:end) = BSal_local(:,BestModSal);
        BSal(tr,1) = FitInfo_Sal.Intercept(BestModSal);
        
        
        %         % spectral median predicting Y
        %         [BSpecMed, FitInfo_SpecMed]=lassoglm(XSpecMed,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
        %         % find the model with the minimum of deviance (best lambda)
        %         [BestDevSpecMed(tr,dd),BestModSpecMed] = min(FitInfo_SpecMean.Deviance);
        
        % amplitude predicting Y, is a null model for the former 2 models
        [BAmp_local, FitInfo_Amp]=lassoglm(XAmp,Y,ParamModel.DISTR,'Alpha', ParamModel.Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',1,'LambdaRatio',ParamModel.LAMBDARATIO);
        % find the model with the minimum of deviance (best lambda)
        [BestDevAmp(tr),BestModAmp] = min(FitInfo_Amp.Deviance);
        
        BAmp(tr,2:end) = BAmp_local(:,BestModAmp);
        BAmp(tr,1) = FitInfo_Amp.Intercept(BestModAmp);
        
        % null model
        MDL_null=fitglm(ones(size(XAmp,1),1),Y,'Distribution',ParamModel.DISTR,'Link',ParamModel.LINK,'Intercept',false);
        % find the model with the minimum of deviance (best lambda)
        %     [BestDevNull(tr),BestModNull] = min(FitInfo_Amp.Deviance);
        %
        %     BNull(tr,2:end) = Bnull_local(:,BestModNull);
        %     BNull(tr,1) = FitInfo_null.Intercept(BestModNull);
        
        BestDevNull(tr) = MDL_null.Deviance;
        BNull(tr) = MDL_null.Coefficients.Estimate;
        
        if DatFig
            TimeBinsX = -Delay : (Win-Delay);
            figure()
            plot(TimeBinsX,BSpecMean(:,BestModSpecMean),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Spectral Mean Dev = %.1f',BestDevSpecMean(tr)))
            hold on
            plot(TimeBinsX,BSal(:,BestModSal),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Saliency Dev = %.1f',BestDevSal(tr)))
            hold on
            %         plot(TimeBinsX,BSpecMed(:,BestModSpecMed),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Spectral Median Dev = %.1f',BestDevSpecMed(tr,dd)))
            hold on
            plot(TimeBinsX,BAmp(:,BestModAmp),'LineStyle','-', 'LineWidth',2,'DisplayName', sprintf('Amp Dev = %.1f',BestDevAmp(tr)))
            
            
            legend('show')
            % XTick = get(gca, 'XTickLabel');
            % XTick = cellfun(@str2double, XTick) * Win;
            % set(gca,'XTickLabel',XTick)
            xlabel('Time (ms)')
            title(sprintf('Poisson Ridge regression on Acoustic Features Delay = %dms Time resolution=%dms', Delay, TR))
            hold off
            pause(1)
        end
        TicToc(tr) = toc(TimerLoop);
        fprintf(1,'Cell %d/%d Models with Time resolution %d ms (%d/%d) => done in %.2f minutes \n', cc,NCells,TR, tr, length(TRs), TicToc(tr)/60);
    end
    if DatFig
        
        figure()
        ColorCode = get(groot, 'DefaultAxesColorOrder');
        subplot(1,2,1)
        plot(TRs,BestDevSpecMean', 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        hold on
        plot(TRs,BestDevNull, 'Color',ColorCode(4,:),'LineWidth',2, 'LineStyle','--')
        legend({'Amp + SpecMean' 'Amp + Sal' 'Amp' 'Null'})
        legend('Autoupdate','off')
        ylabel('Deviance')
        xlabel('Time resolution in ms')
        YLim = get(gca,'YLim');
        
        subplot(1,2,2)
        plot(TRs,BestDevSpecMean', 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        hold on
        plot(TRs,BestDevNull, 'Color',ColorCode(4,:),'LineWidth',2, 'LineStyle','--')
        ylabel('Deviance')
        xlabel('Time resolution in ms')
        xlim([0 2])
        ylim([0.9*YLim(2) YLim(2)])
        hold off
        
        
        figure()
        subplot(1,2,1)
        plot(TRs, BestDevAmp-BestDevSpecMean , 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevAmp - BestDevSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs,BestDevNull -BestDevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
        ylabel('Difference of Deviance')
        xlabel('Time resolution in ms')
        
        subplot(1,2,2)
        plot(TRs, (BestDevAmp-BestDevSpecMean)./BestDevAmp , 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(TRs,(BestDevAmp - BestDevSal)./BestDevAmp,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(TRs,(BestDevNull -BestDevAmp)./BestDevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
        ylabel('R2')
        xlabel('Time resolution in ms')
        
        % Plotting betas
        figure()
        subplot(2,3,1:3)
        bar(TRs,exp([BNull BAmp(:,1) BSpecMean(:,1) BSal(:,1)])*1000)
        xlabel('Time resolution in ms')
        ylabel('exp(intercept) Hz')
        legend({'Null' 'Amp' 'Amp + SpecMean' 'Amp + Sal'},'Location','southoutside','NumColumns',4)
        subplot(2,3,4)
        imagesc(BAmp(:,2:end))
        set(gca,'YTick',1:length(TRs),'YTickLabel',TRs)
        ylabel('Time resolution in ms')
        xlabel('Time alligned to predicted Y (ms)')
        set(gca,'XTick',1:50:Win,'XTickLabel',-Delay:50:(Win-Delay-1))
        title('Amplitude Model')
        colorbar()
        subplot(2,3,5)
        imagesc(BSpecMean(:,2:end))
        set(gca,'YTick',1:length(TRs),'YTickLabel',TRs)
        ylabel('Time resolution in ms')
        xlabel('Time alligned to predicted Y (ms)')
        set(gca,'XTick',1:100:2*Win,'XTickLabel',[-Delay:100:(Win-Delay-1) -Delay:100:(Win-Delay-1)])
        title('Amplitude + SpecMean')
        colorbar()
        subplot(2,3,6)
        imagesc(BSal(:,2:end))
        set(gca,'YTick',1:length(TRs),'YTickLabel',TRs)
        ylabel('Time resolution in ms')
        xlabel('Time alligned to predicted Y (ms)')
        set(gca,'XTick',1:100:2*Win,'XTickLabel',[-Delay:100:(Win-Delay-1) -Delay:100:(Win-Delay-1)])
        title('Amplitude + Saliency')
        colorbar()
    end
    
    
    MotorModels{cc}.TRs = TRs;
    MotorModels{cc}.DevAmp = BestDevAmp;
    MotorModels{cc}.DevAmpSpecMean = BestDevSpecMean;
    MotorModels{cc}.DevAmpSal = BestDevSal;
    MotorModels{cc}.DevNull = BestDevNull;
    MotorModels{cc}.BAmp = BAmp;
    MotorModels{cc}.BAmpSpecMean = BSpecMean;
    MotorModels{cc}.BAmpSal = BSal;
    MotorModels{cc}.BNull = BNull;
    MotorModels{cc}.TicToc = TicToc;
    MotorModels{cc}.MeanY = MeanY;
        
end
%%
save(fullfile(Path,'MotorModelsAllCells'),'MotorModels','CellsPath');
%% Plot issues with no convergence
load(fullfile(Path,'MotorModelsAllCells'),'MotorModels','CellsPath');
TicToc = nan(NCells*11,1);
MeanY = nan(NCells*11,1);
TRs = nan(NCells*11,1);
BestDevAmp = nan(NCells*11,1);
BestDevSpecMean = nan(NCells*11,1);
BestDevSal = nan(NCells*11,1);
BestDevNull = nan(NCells*11,1);

F10=figure(10);
F20 = figure(20);
ColorCode = get(groot, 'DefaultAxesColorOrder');
ii=0;
for cc=1:NCells
    if ~isempty(MotorModels{cc})
        for tt=1:length(MotorModels{cc}.TicToc)
            ii=ii+1;
            TicToc(ii) = MotorModels{cc}.TicToc(tt);
            MeanY(ii) = MotorModels{cc}.MeanY(tt);
            TRs(ii) = MotorModels{cc}.TRs(tt);
            BestDevAmp(ii) = MotorModels{cc}.DevAmp(tt);
            BestDevSpecMean(ii) = MotorModels{cc}.DevAmpSpecMean(tt);
            BestDevSal(ii) = MotorModels{cc}.DevAmpSal(tt);
            BestDevNull(ii) = MotorModels{cc}.DevNull(tt);
        end
        
        figure(10)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevAmpSpecMean', 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevAmpSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevNull, 'Color',ColorCode(4,:),'LineWidth',2, 'LineStyle','--')
        if cc==1
            legend({'Amp + SpecMean' 'Amp + Sal' 'Amp' 'Null'})
            legend('Autoupdate','off')
            ylabel('Deviance')
            xlabel('Time resolution in ms')
        end
        
        
        figure(20)
        subplot(1,2,1)
        hold on
        plot(MotorModels{cc}.TRs, MotorModels{cc}.DevAmp-MotorModels{cc}.DevAmpSpecMean , 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevAmp - MotorModels{cc}.DevAmpSal,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,MotorModels{cc}.DevNull -MotorModels{cc}.DevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        if cc==1
            legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
            legend('Autoupdate','off')
            ylabel('Difference of Deviance')
            xlabel('Time resolution in ms')
        end
        
        subplot(1,2,2)
        hold on
        plot(MotorModels{cc}.TRs, (MotorModels{cc}.DevAmp-MotorModels{cc}.DevAmpSpecMean)./MotorModels{cc}.DevAmp , 'Color',ColorCode(1,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,(MotorModels{cc}.DevAmp - MotorModels{cc}.DevAmpSal)./MotorModels{cc}.DevAmp,'Color', ColorCode(2,:), 'LineWidth',2)
        hold on
        plot(MotorModels{cc}.TRs,(MotorModels{cc}.DevNull -MotorModels{cc}.DevAmp)./MotorModels{cc}.DevAmp, 'Color', ColorCode(3,:),'LineWidth',2)
        if cc==1
            legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
            legend('Autoupdate','off')
            ylabel('R2')
            xlabel('Time resolution in ms')
        end
        
    end
end

figure()
subplot(1,2,1)
scatter(MeanY,TicToc/(60*60), 'MarkerFaceColor','k')
xlabel('MeanRate')
ylabel('Model Time consumption (h)')
subplot(1,2,2)
scatter(TRs,TicToc/(60*60), 'MarkerFaceColor','k')
xlabel('Models Time Resolution (ms)')
ylabel('Model Time consumption (h)')


figure()
ColorCode = get(groot, 'DefaultAxesColorOrder');
%subplot(1,2,1)
scatter(TRs,BestDevSpecMean, 'MarkerFaceColor',ColorCode(1,:), 'MarkerEdgeColor',ColorCode(1,:))
hold on
scatter(TRs,BestDevSal,'MarkerFaceColor',ColorCode(2,:), 'MarkerEdgeColor',ColorCode(2,:))
hold on
scatter(TRs,BestDevAmp, 'MarkerFaceColor',ColorCode(3,:), 'MarkerEdgeColor',ColorCode(3,:))
hold on
scatter(TRs,BestDevNull, 'MarkerFaceColor',ColorCode(4,:), 'MarkerEdgeColor',ColorCode(4,:))
legend({'Amp + SpecMean' 'Amp + Sal' 'Amp' 'Null'})
legend('Autoupdate','off')
ylabel('Deviance')
xlabel('Time resolution in ms')




figure()
subplot(1,2,1)
scatter(TRs, BestDevAmp-BestDevSpecMean , 'MarkerFaceColor',ColorCode(1,:), 'MarkerEdgeColor',ColorCode(1,:))
hold on
scatter(TRs,BestDevAmp - BestDevSal,'MarkerFaceColor',ColorCode(2,:), 'MarkerEdgeColor',ColorCode(2,:))
hold on
scatter(TRs,BestDevNull -BestDevAmp, 'MarkerFaceColor',ColorCode(3,:), 'MarkerEdgeColor',ColorCode(3,:))
legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
ylabel('Difference of Deviance')
xlabel('Time resolution in ms')

subplot(1,2,2)
scatter(TRs, (BestDevAmp-BestDevSpecMean)./BestDevAmp , 'MarkerFaceColor',ColorCode(1,:), 'MarkerEdgeColor',ColorCode(1,:))
hold on
scatter(TRs,(BestDevAmp - BestDevSal)./BestDevAmp,'MarkerFaceColor',ColorCode(2,:), 'MarkerEdgeColor',ColorCode(2,:))
hold on
scatter(TRs,(BestDevNull -BestDevAmp)./BestDevAmp,'MarkerFaceColor',ColorCode(3,:), 'MarkerEdgeColor',ColorCode(3,:))
legend({'SpecMean contribution' 'Sal contribution' 'Amp contribution'})
ylabel('R2')
xlabel('Time resolution in ms')
%%
        





% %% Old plots of deviances
% CLIM = nan(4,2);
% figure()
% imagesc(Delay,TRs,BestDevSal)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance pitch saliency model')
% colorbar()
% CLIM(2,:) = get(gca, 'clim');
% 
% % figure()
% % imagesc(Delays,TRs,BestDevSpecMed)
% % xlabel('Delays to voc onset in ms')
% % ylabel('Time resolution in ms')
% % title('Deviance spectral median model')
% 
% figure()
% imagesc(Delay,TRs,BestDevAmp)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance Amplitude model')
% colorbar()
% CLIM(3,:) = get(gca, 'clim');
% 
% figure()
% imagesc(Delay,TRs,BestDevNull)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance Null model')
% colorbar()
% CLIM(4,:) = get(gca, 'clim');
% 
% 
% figure()
% subplot(2,3,1)
% imagesc(BestDevSpecMean)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance spectral mean model')
% colorbar()
% %caxis([0 1000])
% caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% 
% subplot(2,3,2)
% imagesc(BestDevSal)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance pitch saliency model')
% colorbar()
% caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% %caxis([0 1000])
% 
% % figure()
% % imagesc(Delays,TRs,BestDevSpecMed)
% % xlabel('Delays to voc onset in ms')
% % ylabel('Time resolution in ms')
% % title('Deviance spectral median model')
% 
% subplot(2,3,3)
% imagesc(BestDevAmp)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance Amplitude model')
% colorbar()
% caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% %caxis([0 1000])
% 
% subplot(2,3,5)
% imagesc(BestDevNull)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Deviance Null model')
% colorbar()
% caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% 
% figure()
% subplot(1,3,1)
% imagesc(BestDevAmp - BestDevSpecMean)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Goodness of fit: Deviance Amp - SpecMean')
% colorbar()
% % caxis([0 1000])
% % %caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% 
% subplot(1,3,2)
% imagesc(BestDevAmp - BestDevSal)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Goodness of fit: Deviance Amp - pitch saliency')
% colorbar()
% % %caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% % caxis([0 1000])
% 
% % figure()
% % imagesc(Delays,TRs,BestDevSpecMed)
% % xlabel('Delays to voc onset in ms')
% % ylabel('Time resolution in ms')
% % title('Deviance spectral median model')
% 
% subplot(1,3,3)
% imagesc(BestDevNull - BestDevAmp)
% set(gca, 'XTick',1:length(Delay),'XTickLabel',Delay)
% set(gca, 'YTick',1:length(TRs),'YTickLabel',TRs)
% xlabel('Delays to voc onset in ms')
% ylabel('Time resolution in ms')
% title('Goodness of fit: Deviance Null -  Amplitude')
% colorbar()
% %caxis([min(CLIM(:,1)) max(CLIM(:,2))])
% caxis([0 1000])

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

function [XPerStim] = get_x(BioSound, Duration, Win,Delay,Feature)
% calculate the cell array of acoustic data for each vocalization. For each
% voc, XPerStim is a matrix where each column correspond to the values of
% the acoustic feature in the window Win starting Delay ms before time t of
% the neural response Y. Same as the neural response, the window slides by
% a pace of 1ms from one column of the matrix to the next, starting
% Win-Delay before the vocalization onset.

% Matrix of the acoustic features
XPerStim = cell(1,length(Duration));
for stim = 1:length(Duration)
    % Time points of the corresponding neural response
    TimeBinsY = -(Win-Delay) : 1 : (Delay + Duration(stim));
    
    % Get ready the stim acoustic features that was sampled at 1000Hz
    FeatureVal = BioSound{stim}.(sprintf('%s',Feature));
    %         FeatureValResampled = resample(FeatureVal,1,TR); % Win in ms so resampling from 1 value per ms to 1/Win value per ms or 1 value per win ms
    %         PadFeature4Conv = nanmean(FeatureVal)*ones(1,(length(Expwav)-1)/2);
    %         FeatureValResampled = conv([PadFeature4Conv FeatureVal PadFeature4Conv], Expwav, 'valid');
    ZPaddFeature = nan(1,Win + Duration(stim) + Win);
    ZPaddFeature(Win+(1:length(FeatureVal)))=FeatureVal;
    
    XPerStim{stim} = zeros(Win,length(TimeBinsY)-1);
    for tt=1:(length(TimeBinsY)-1)
        XPerStim{stim}(:,tt) = ZPaddFeature((1:Win) + tt-1);
    end
end

end

function [YPerStim] = get_y(SAT, Duration, Win,Delay,TR)
% Calculate the time varying rate applying a gaussian window TR on the
% spike pattern. The spike pattern considered starts Win-Delay ms
% before the onset of the vocalization and stops Delay ms after the
% offset of the vocalization
YPerStim = cell(1,length(Duration));
% Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
nStd =4;
Tau = (TR/2);
T_pts = (0:2*nStd*Tau) -nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to 4 std away on each side
Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
Expwav = Expwav./sum(Expwav);

% Loop through the stimuli and fill in the matrix
for stim=1:length(Duration)
    % Time slots for the neural response
    TimeBinsY = -(Win-Delay) : (Delay + Duration(stim));
    SpikePattern = zeros(1,length(TimeBinsY)-1);
    for isp = 1:length(SAT{stim})
        SpikeInd = round(SAT{stim}(isp));
        if (SpikeInd>=-(Win-Delay)) && (SpikeInd<(Delay + Duration(stim)))
            SpikePattern(SpikeInd + (Win-Delay) +1) = SpikePattern(SpikeInd + (Win-Delay) +1) +1;
        end
    end
    YPerStim{stim} = conv(SpikePattern, Expwav,'same');
    %         % Time slots for the neural response
    %         TimeBinsY = -(Delay) : TR: (Delay + Duration(stim));
    %         YPerStim{stim} = nan(1,length(TimeBinsY)-1);
    %         for tt=1:(length(TimeBinsY)-1)
    %             % Find the number of spikes
    %             YPerStim{stim}(tt) = sum( (SAT{stim}>=TimeBinsY(tt)) .* (SAT{stim}<TimeBinsY(tt+1)));
    %         end
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

function plotxyfeatures(BioSound,YPerStim,XPerStim,TR,Win,Delay,Duration,F_high,FeatureName)
if nargin<6
    FeatureName = 'Acoustic Feature';
end
DBNOISE =60;
f_low = 0;
close all
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
    StimFeature = [XPerStim{stim}(:,1)' XPerStim{stim}(end,2:end)];
    MaxT = floor(max(double(BioSound{stim}.to)*1000));
    XStimFeature = [-Win:-1 double(BioSound{stim}.to)*1000 MaxT+(1:(Win))];
    plot(XStimFeature(1:length(StimFeature)), StimFeature, 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    %     IndMax = floor(max(double(BioSound{stim}.to)*1000)/TR);
    %     plot(TR/2 + (0:TR:((IndMax-1)*TR)), StimFeature(1:IndMax), 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    ylabel(sprintf('%s', FeatureName))
    if strcmp(FeatureName,'Spectral Mean')
        ylim(v_axis(3:4))
    end
    xlim([-Win Win+Duration(stim)])
    hold off
    
    subplot(2,1,2)
    yyaxis left
    bar(-(Win-Delay) : (Delay + Duration(stim)-1),YPerStim{stim})
    %     bar(TR/2+(-(Delay) : TR : (Delay + Duration(stim)-TR)),YPerStim{stim})
    ylabel(sprintf('Spikes/ms (%d ms Gauss win)', TR))
    ylim([0 1])
    xlabel('Time (ms)')
    hold on
    line([-(Win-Delay) Duration(stim)+Delay], [0.975 0.975], 'Color',[0 0.4470 0.7410], 'LineWidth',12)
    hold on
    text(0,0.975,'Neural response', 'Color', [1 1 1])
    hold on
    yyaxis right
    plot(XStimFeature(1:length(StimFeature)), StimFeature, 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    %     plot(TR/2 + (0:TR:((IndMax-1)*TR)), StimFeature(1:IndMax), 'Color',ColorCode(2,:),'LineStyle','-', 'LineWidth',2)
    ylabel(sprintf('%s', FeatureName))
    xlim([-Win Win+Duration(stim)])
    hold off
    pause(1)
end
end



