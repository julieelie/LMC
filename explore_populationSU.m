% Path 2 Data
Path2AllData = '/Volumes/server_home/users/JulieE/LMC_HoHa/logger/';
addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'))
ADBitVolts_sorting=500/32767; %Conversion value used in mat2ntt for proper ploting under spikesort3D in AD count
% Let's find the date of interest in the Google sheet
Path2RecordingTable = '/Users/elie/Google Drive/BatmanData/RecordingLogs/recording_logs.xlsx';
[~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:AF200','basic');
Header = RecTableData(1,:);
ND_IDCol = find(contains(Header, 'Neural data voc in operant + PSTH'));
BatIDCol = find(contains(Header, 'Bat'));
NLCol = find(contains(Header, 'NL'));
Dates_Ind = find(cell2mat(RecTableData(2:end,ND_IDCol))==1)+1;
Dates = cell2mat(RecTableData(Dates_Ind,1));
ThreshSig=50;
%% Loop through dates and gather data
% For each day, each logger, each unit find the average spike rate
PLOT=1;
MeanSR4Calulations = 1; % Values in Hz of the mean SR for doing calculations
NDays = length(Dates);
HypNSU = 200;
TimeStep = 10*60; % Time resolution in seconds at which the spike rate shoule be calculated
MaxRecTime = 60*60*6;
TimePoints = 0:TimeStep:MaxRecTime;
SpikeRate = nan(HypNSU,length(TimePoints)-1);% will contain spike rate in Hz
KDE = nan(HypNSU,length(TimePoints)-1);% will contain spike rate in Hz
KDE_error = cell(HypNSU,1);% will contain spike rate error in Hz
Peak2Peak = nan(HypNSU,length(TimePoints)-1);% will contain mean spike amplitude along time
SpikeLife = nan(HypNSU,1); % Duration in minutes of period containing a spike from the cell
SD_SR = nan(HypNSU,1);
Mean_SR = nan(HypNSU,1); % average mean SR in Hz of each cell only counting time windows of 10min where there is at least one spike (cell still active)
Bat_ID = nan(HypNSU,1);
Date_ID = nan(HypNSU,1);
Tetrode_ID = nan(HypNSU,1);
SS_ID = cell(HypNSU,1);
Mean_Peak2Peak = nan(HypNSU,1);
SpikeShape = cell(HypNSU,4);
Alpha = 0.01; % Detection thershold for significant spike rate modulation before correction for false positive
SignWinVoc = nan(HypNSU,1); % number of time bin in the KDE when producing the vocalizations above the mean baseline firing rate by twice the SD of the baseline firing rate
SignWinHear = nan(HypNSU,1);% number of time bin in the KDE when hearing the vocalizations above the mean baseline firing rate by twice the SD of the baseline firing rate
ProdCellGood = nan(HypNSU,1);

CellCount=0;
for dd=1:NDays
    Date = Dates(dd);
    fprintf(1,'%d\n',Date)
    Path2Data = fullfile(Path2AllData, ['20' num2str(Date)]);
    % Get the number of loggers
    Logger_dirs = dir(fullfile(Path2Data, '*ogger*'));
    Logger_dirs=Logger_dirs([Logger_dirs.isdir]);
    NLogger = length(Logger_dirs);
    % Identify the type of logger and extract neural data
    LoggerType = cell(NLogger,1);
    % LOAD the data of KDE to do a test for vocalization production
    FileData = dir(fullfile(Path2Data, sprintf('%d_*_VocExtractData_*', Date)));
    load(fullfile(FileData(1).folder,FileData(1).name), 'SpikeTimesVoc')
    Ind_ = strfind(FileData(1).name, '_');
    Delay = str2double(FileData(1).name(Ind_(end) + (1:3)));
    for ll=1:NLogger
        fprintf(1,'%s\n',Logger_dirs(ll).name)
        LData_folder = fullfile(Logger_dirs(ll).folder, Logger_dirs(ll).name,'extracted_data');
        LDir = dir(fullfile(LData_folder, '*CSC*.mat'));
        LData = load(fullfile(LDir(1).folder, LDir(1).name), 'logger_type', 'logger_serial_number');
        LoggerType{ll}  = LData.logger_type;
        if strcmp(LoggerType{ll}, 'Mous') || strcmp(LoggerType{ll}, 'Rat')
            SU_files = dir(fullfile(LData_folder, '*TT*SS*.mat'));
            NSU = length(SU_files);
            %Get the Bat-ID for that logger
            RowData = find((cell2mat(RecTableData(2:end,1))== Date)) +1;
            DataInfo = RecTableData(RowData,:);
            NL_ID = cell2mat(DataInfo(NLCol));
            NLCol_local = NLCol(NL_ID==str2double(LData.logger_serial_number));
            Bat_ID_local = DataInfo{BatIDCol(find(BatIDCol<NLCol_local,1,'last'))};
            if Bat_ID_local==11689 % Only run on Hodor

                for uu=1:NSU
                    CellCount = CellCount+1;
                    fprintf(1,'Cell %d\n',CellCount)
                    Cell = load(fullfile(SU_files(uu).folder, SU_files(uu).name)); % This load the Spike_arrival_times in micro seconds and the Spike_snippets in microVolts, ceiled at 500uV for spike sorting projection
                    SpikeTimes = (Cell.Spike_arrival_times-Cell.Spike_arrival_times(1))*10^-6; % Spike arrival times centered to the first spike, in seconds

                    % Find the best channel to calculate the peak2peak value
                    % (channel with largest spike trace)
                    P2P_all = nan(size(Cell.Spike_snippets,3),4);
                    for cc=1:4 % calculate for each channel the peak2peak for all spikes
                        P2P_all(:,cc)=reshape(max(Cell.Spike_snippets(:,cc,:),[],1) - min(Cell.Spike_snippets(:,cc,:),[],1),size(Cell.Spike_snippets,3),1);
                    end
                    [~,I]=max(P2P_all,[],2);
                    Best_c = nan(1,4);
                    for cc=1:4
                        Best_c(cc) = sum(I==cc);
                    end
                    [~,Best_c] = max(Best_c); % This is the channel with the largest spike


                    for tt=1:(length(TimePoints)-1)
                        Spike_local = (SpikeTimes>=TimePoints(tt)) .* (SpikeTimes<=TimePoints(tt+1));
                        SpikeRate(CellCount,tt) = sum(Spike_local)/TimeStep; % this is the local spike rate in Hertz
                        SpikeInd = find(Spike_local);
                        % calculate for each channel the average peak2peak for the spikes occuring during that time slot
                        Peak2Peak(CellCount,tt) = ADBitVolts_sorting.*mean(max(Cell.Spike_snippets(:,Best_c,SpikeInd),[],1) - min(Cell.Spike_snippets(:,Best_c,SpikeInd),[],1));
                    end



                    % calculate average statistics for that cell
                    SpikeLife(CellCount) = (TimePoints(find(~isnan(Peak2Peak(CellCount,:)),1,'Last'))- TimePoints(find(~isnan(Peak2Peak(CellCount,:)),1,'First')))/60;
                    [KDE(CellCount,:),~,KDE_error{CellCount}] = kde_wrapper(SpikeTimes,TimePoints(2:end)-TimeStep/2,1/TimeStep);
                    Mean_Peak2Peak(CellCount) = mean(P2P_all(:,Best_c));
                    Mean_SR(CellCount) = mean(SpikeRate(CellCount,find(~isnan(Peak2Peak(CellCount,:)))));
                    SD_SR(CellCount) = std(SpikeRate(CellCount,find(~isnan(Peak2Peak(CellCount,:)))));
                    Bat_ID(CellCount) = Bat_ID_local;
                    Date_ID(CellCount) = str2double(['20' num2str(Date)]);
                    TTInd = strfind(SU_files(uu).name, 'TT');
                    Tetrode_ID(CellCount) = str2double(SU_files(uu).name(TTInd+2));
                    SSInd = strfind(SU_files(uu).name, 'SS');
                    SS_ID{CellCount} = SU_files(uu).name(SSInd+(3:4));

                    % extract the average spike shape for that cell
                    for cc=1:4
                        SpikeShape{CellCount,cc}=nan(2,size(Cell.Spike_snippets,1));
                        SpikeShape{CellCount,cc}(1,:) = mean(Cell.Spike_snippets(:,cc,:),3).*ADBitVolts_sorting;
                        SpikeShape{CellCount,cc}(2,:) = std(Cell.Spike_snippets(:,cc,:),0,3).*ADBitVolts_sorting;
                    end
                    
                    
                    % Test if the cell has a modulated spike rate compared to
                    % averaged before vocalizations (one time window with KDE
                    % away from the mean by twice SD). Rate is in spike per ms
                    % here. Test only if the average spike rate is  above 0.1Hz (1 spike every 10 seconds);
                    ST = SpikeTimesVoc.(sprintf('Logger%s',LData.logger_serial_number));
                    LogicalInd=logical((ST.Sum_Psth_KDEfiltered_VocCall{uu,1}>-Delay) .* (ST.Sum_Psth_KDEfiltered_VocCall{uu,1}<mean(ST.VocDuration+Delay)));
                    Fig0=figure(1);
                    clf(Fig0)
                    plot(ST.Sum_Psth_KDEfiltered_VocCall{uu,1}(LogicalInd), ST.Sum_Psth_KDEfiltered_VocCall{uu,2}(LogicalInd), '-k', 'LineWidth',2)
                    hold on
                    plot([min(ST.Sum_Psth_KDEfiltered_VocCall{uu,1}(LogicalInd)) max(ST.Sum_Psth_KDEfiltered_VocCall{uu,1}(LogicalInd))], MeanBaselineVoc/1000.* ones(2,1), ':r', 'LineWidth',2)
                    ylabel('Spike Rate (Hz)')
                    xlabel('Time (ms)')
                    legend('Rate','Background')
                    title(sprintf('Test for modulated spike rate M%d %d TT%d SS%s SU%d', Bat_ID(CellCount),Date_ID(CellCount), Tetrode_ID(CellCount), SS_ID{CellCount},uu))
                    
                    hold off
                    if Mean_SR(CellCount)>MeanSR4Calulations
                        if ~isempty(ST.Sum_Psth_KDEfiltered_VocBaseline) && iscell(ST.Sum_Psth_KDEfiltered_VocBaseline)
                            MeanBaselineVoc = mean(ST.Sum_Psth_KDEfiltered_VocBaseline{uu,2})*1000;% Average baseline rate in Hz
                            if MeanBaselineVoc>MeanSR4Calulations
                                Zscore = abs(ST.Sum_Psth_KDEfiltered_VocCall{uu,2}(LogicalInd) - MeanBaselineVoc/1000)/std(ST.Sum_Psth_KDEfiltered_VocBaseline{uu,2});
                                P_value = normcdf(Zscore, 'upper');
                                PV_sorted = sort(P_value);
                                
                                
                                % False Rate detection correction,
                                % Benjamini-Hochberg procedure
                                SignWinVoc(CellCount)=sum(PV_sorted<(Alpha/length(PV_sorted).*(1:length(PV_sorted))));
                            end
                        end
                        if ~isempty(ST.Sum_Psth_KDEfiltered_HearBaseline) && iscell(ST.Sum_Psth_KDEfiltered_HearBaseline)
                            MeanBaselineHear = mean(ST.Sum_Psth_KDEfiltered_HearBaseline{uu,2})*1000;% Average baseline rate in Hz
                            if MeanBaselineHear>MeanSR4Calulations
                                LogicalInd=logical((ST.Sum_Psth_KDEfiltered_HearCall{uu,1}>-Delay) .* (ST.Sum_Psth_KDEfiltered_HearCall{uu,1}<mean(ST.HearDuration(ST.HearOnlyInd)+Delay)));
                                Zscore = abs(ST.Sum_Psth_KDEfiltered_HearCall{uu,2}(LogicalInd) - MeanBaselineHear/1000)/std(ST.Sum_Psth_KDEfiltered_HearBaseline{uu,2});
                                P_value = normcdf(Zscore, 'upper');
                                PV_sorted = sort(P_value);
                                % False Rate detection correction,
                                % Benjamini-Hochberg procedure
                                SignWinHear(CellCount)=sum(PV_sorted<(Alpha/length(PV_sorted).*(1:length(PV_sorted))));
                            end
                        end

                    end

                    if PLOT
                        % Plot the average spike rate over time, the average spike
                        % amplitude over time and the mean snippets
                        Fig1=figure(2);
                        clf(Fig1)
                        yyaxis left
                        plot((TimePoints(2:end)-TimeStep/2)/60,SpikeRate(CellCount,:),'b-', 'LineWidth',2)
                        hold on
                        shadedErrorBar((TimePoints(2:end)-TimeStep/2)/60,KDE(CellCount,:),KDE_error{CellCount},{'b--', 'LineWidth',2})
                        ylabel('Spike rate (Hz)')
                        xlabel('Time (min)')
                        ylim([0 max(MeanSR4Calulations/10,Fig1.Children.YLim(2))])
                        hold on
                        if Fig1.Children.YLim(2)>MeanSR4Calulations
                            hline(MeanSR4Calulations, 'k:', 'Spike rate threshold for calculating significance')
                        end
                        yyaxis right
                        plot((TimePoints(2:end)-TimeStep/2)/60,Peak2Peak(CellCount,:),'r-', 'LineWidth',2)
                        ylabel('Spike Amplitude (uV)')
                        Fig1.Children.YLim(1)=0;
                        title(sprintf('M%d %d TT%d SS%s SU%d SignifVoc=%d SignifHear=%d', Bat_ID(CellCount),Date_ID(CellCount), Tetrode_ID(CellCount), SS_ID{CellCount},uu,SignWinVoc(CellCount),SignWinHear(CellCount)))
                        hold off
                        Fig2=figure(3);
                        clf(Fig2)
                        Fig2YLim = nan(2,4);
                        for cc=1:4
                            subplot(2,2,cc)
                            shadedErrorBar([],SpikeShape{CellCount,cc}(1,:), SpikeShape{CellCount,cc}(2,:), {'Color','k','LineWidth',2})
                            ylabel('Voltage (uVolt)')
                        end
                        for cc=1:4
                            Fig2YLim(:,cc) = Fig2.Children(cc).YLim;
                        end
                        Fig2YLim_final(1) = min(Fig2YLim(1,:));
                        Fig2YLim_final(2) = max(Fig2YLim(2,:));
                        for cc=1:4
                            Fig2.Children(cc).YLim = Fig2YLim_final;
                        end
                        sgtitle(sprintf('Spike shape M%d %d TT%d SS%s SU%d', Bat_ID(CellCount),Date_ID(CellCount), Tetrode_ID(CellCount), SS_ID{CellCount},uu))
                        print(Fig1,fullfile(Path2Data,sprintf('%d_Logger%s_TT%d_SS%s_SU%d_Rate.pdf', Date, LData.logger_serial_number,Tetrode_ID(CellCount),SS_ID{CellCount},uu)),'-dpdf')
                        print(Fig2,fullfile(Path2Data,sprintf('%d_Logger%s_TT%d_SS%s_SU%d_Snippets.pdf', Date, LData.logger_serial_number,Tetrode_ID(CellCount),SS_ID{CellCount},uu)),'-dpdf')
    %                     if SignWinVoc(CellCount)>0 || SignWinHear(CellCount)>0
        %                     pause()
        %                 end
                        pause(5)
                    end
                    if SignWinVoc(CellCount)>ThreshSig
                        ProdCellGood(CellCount)=input('Good Production Cell?');
                    end
                end
            end
            
        end
    end
end


% Plot the histogram of cells Mean_SR
figure()
[~,edges] = histcounts(log10(Mean_SR));
histogram(Mean_SR,10.^edges)
set(gca, 'XScale','log')
ylabel('# Cells')
xlabel('Average rate (Hz)')

% Scatter plot of mean_SR and Mean amplitude of the spikes

figure()
subplot(2,2,1)
plot(Mean_SR, Mean_Peak2Peak, 'k.','MarkerSize',10)
hold on
plot(Mean_SR(SignWinVoc>ThreshSig), Mean_Peak2Peak(SignWinVoc>ThreshSig), 'r.','MarkerSize',10)
set(gca, 'XScale','log')
ylabel('Spike Amplitude uV')
xlabel('Average rate (Hz)')
hold off

subplot(2,2,2)
plot(Mean_SR, Mean_Peak2Peak, 'k.','MarkerSize',10)
hold on
plot(Mean_SR(SignWinHear>ThreshSig), Mean_Peak2Peak(SignWinHear>ThreshSig), 'c.','MarkerSize',10)
set(gca, 'XScale','log')
ylabel('Spike Amplitude uV')
xlabel('Average rate (Hz)')
hold off

% Scatter plot of mean_SR and Recording duration
subplot(2,2,3)
plot(Mean_SR, SpikeLife, 'k.','MarkerSize',10)
hold on
plot(Mean_SR(SignWinVoc>ThreshSig), SpikeLife(SignWinVoc>ThreshSig), 'r.','MarkerSize',10)
set(gca, 'XScale','log')
ylabel('Recording duration (min)')
ylim([0 400])
xlabel('Average rate (Hz)')
hold off

subplot(2,2,4)
plot(Mean_SR, SpikeLife, 'k.','MarkerSize',10)
hold on
plot(Mean_SR(SignWinHear>ThreshSig), SpikeLife(SignWinHear>ThreshSig), 'c.','MarkerSize',10)
set(gca, 'XScale','log')
ylim([0 400])
ylabel('Recording duration (min)')
xlabel('Average rate (Hz)')
hold off
save('/Volumes/server_home/users/JulieE/LMC_HoHa/PopulationAnalysis.mat')


%% Look at cell one by one
load('/Volumes/server_home/users/JulieE/LMC_HoHa/PopulationAnalysis.mat')
% First voc prod modulated
% ProdCell = find(SignWinVoc>ThreshSig);
% for cc=1:length(ProdCell)
%     fprintf(1,'Date: %d\n',Date_ID(ProdCell(cc)))
%     fprintf(1,'Bat: %d\n',Bat_ID(ProdCell(cc)))
%     fprintf(1,'T%d\n',Tetrode_ID(ProdCell(cc)));
%     fprintf(1,'SS%s\n',SS_ID{ProdCell(cc)});
%     % Scatter plot of mean_SR and Mean amplitude of the spikes
%     figure()
%     subplot(2,1,1)
%     plot(Mean_SR, Mean_Peak2Peak, 'k.','MarkerSize',10)
%     hold on
%     plot(Mean_SR(ProdCell(cc)), Mean_Peak2Peak(ProdCell(cc)), 'r.','MarkerSize',10)
%     set(gca, 'XScale','log')
%     ylabel('Spike Amplitude uV')
%     xlabel('Average rate (Hz)')
%     hold off
%     
%     % Scatter plot of mean_SR and Recording duration
%     subplot(2,1,2)
%     plot(Mean_SR, SpikeLife, 'k.','MarkerSize',10)
%     hold on
%     plot(Mean_SR(ProdCell(cc)), SpikeLife(ProdCell(cc)), 'r.','MarkerSize',10)
%     set(gca, 'XScale','log')
%     ylabel('Recording duration (min)')
%     ylim([0 400])
%     xlabel('Average rate (Hz)')
%     hold off  
% end
ProdCellGoodInd = find(ProdCellGood==1);
% Good Prod Cell: 121,122,123, 103,102,99,58

% Then Hear only modulated cell
HearCell = find((SignWinHear>ThreshSig).*~(SignWinVoc>ThreshSig));
HearCellGood = nan(length(HearCell),1);
for cc=1:length(HearCell)
    fprintf(1,'\n\n')
    fprintf(1,'Date: %d\n',Date_ID(HearCell(cc)))
    fprintf(1,'Bat: %d\n',Bat_ID(HearCell(cc)))
    fprintf(1,'T%d\n',Tetrode_ID(HearCell(cc)));
    fprintf(1,'SS%s\n',SS_ID{HearCell(cc)});
    % Scatter plot of mean_SR and Mean amplitude of the spikes
    figure()
    subplot(2,1,1)
    plot(Mean_SR, Mean_Peak2Peak, 'k.','MarkerSize',10)
    hold on
    plot(Mean_SR(HearCell(cc)), Mean_Peak2Peak(HearCell(cc)), 'c.','MarkerSize',10)
    set(gca, 'XScale','log')
    ylabel('Spike Amplitude uV')
    xlabel('Average rate (Hz)')
    hold off
    
    % Scatter plot of mean_SR and Recording duration
    subplot(2,1,2)
    plot(Mean_SR, SpikeLife, 'k.','MarkerSize',10)
    hold on
    plot(Mean_SR(HearCell(cc)), SpikeLife(HearCell(cc)), 'c.','MarkerSize',10)
    set(gca, 'XScale','log')
    ylim([0 400])
    ylabel('Recording duration (min)')
    xlabel('Average rate (Hz)')
    hold off
    HearCellGood=input('GoodCell?');
end

HearCellGood = find(HearCellGood);
% Only 68 is clearly modulated by vocalization hearing


% Scatter plot of mean_SR and Mean amplitude of the spikes
figure()
subplot(2,2,1)
plot(Mean_SR, Mean_Peak2Peak, 'k.','MarkerSize',10)
hold on
plot(Mean_SR(ProdCellGoodInd), Mean_Peak2Peak(ProdCellGoodInd), 'r.','MarkerSize',10)
set(gca, 'XScale','log')
ylabel('Spike Amplitude uV')
xlabel('Average rate (Hz)')
hold off

subplot(2,2,2)
plot(Mean_SR, Mean_Peak2Peak, 'k.','MarkerSize',10)
hold on
if ~isempty(HearCellGood)
    plot(Mean_SR(HearCellGood), Mean_Peak2Peak(HearCellGood), 'c.','MarkerSize',10)
end
set(gca, 'XScale','log')
ylabel('Spike Amplitude uV')
xlabel('Average rate (Hz)')
hold off

% Scatter plot of mean_SR and Mean amplitude of the spikes
subplot(2,2,3)
plot(Mean_SR, SpikeLife, 'k.','MarkerSize',10)
hold on
plot(Mean_SR(ProdCellGoodInd), SpikeLife(ProdCellGoodInd), 'r.','MarkerSize',10)
set(gca, 'XScale','log')
ylabel('Recording duration (min)')
ylim([0 400])
xlabel('Average rate (Hz)')
hold off

subplot(2,2,4)
plot(Mean_SR, SpikeLife, 'k.','MarkerSize',10)
hold on
if ~isempty(HearCellGood)
    plot(Mean_SR(HearCellGood), SpikeLife(HearCellGood), 'c.','MarkerSize',10)
end
set(gca, 'XScale','log')
ylim([0 400])
ylabel('Recording duration (min)')
xlabel('Average rate (Hz)')
hold off
