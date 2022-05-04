addpath(genpath('/Users/elie/Documents/CODE/GitHub/operant_bats'))
addpath(genpath('/Users/elie/Documents/CODE/GitHub/GeneralCode'))
addpath(genpath('/Users/elie/Documents/CODE/GitHub/LMC'))
addpath(genpath('/Users/elie/Documents/CODE/GitHub/Kilosort2'))
addpath(genpath('/Users/elie/Documents/CODE/GitHub/LoggerDataProcessing'))
addpath(genpath('/Users/elie/Documents/CODE/GitHub/SoundAnalysisBats'))
Path2RecordingTable = '/Users/elie/Google Drive/Mon Drive/BatmanData/RecordingLogs/recording_logs.xlsx';
BasePath = '/Volumes/JulieE8T';
%% RUN audio data extraction for the operant tests
% BasePath = '/Volumes/server_home/users/JulieE/LMC';
ListOfPaths = gather_operant_datapath(BasePath);


Path2Run = find(contains(ListOfPaths, 'HoHa'));
% Path2Run(contains(ListOfPaths(Path2Run), '20190207'))=[]; % No Neural Data in Hodor
% Path2Run(contains(ListOfPaths(Path2Run), '20190213'))=[]; % No Neural Data in Hodor
% Path2Run(contains(ListOfPaths(Path2Run), '20190214'))=[]; % No Neural Data in Hodor

%  Path2Run = find(contains(ListOfPaths, 'CoEd'));
%  Path2Run(contains(ListOfPaths(Path2Run), '20190703'))=[];
%  Path2Run(contains(ListOfPaths(Path2Run), '20190709'))=[];
%  Path2Run(contains(ListOfPaths(Path2Run), '190701_0951'))=[];
%  Path2Run(contains(ListOfPaths(Path2Run), '2020'))=[];
%  Path2Run(contains(ListOfPaths(Path2Run), '201905'))=[];% No neural data
%  Path2Run(contains(ListOfPaths(Path2Run), '190605_1406'))=[]; % No vocalization
%  Path2Run(contains(ListOfPaths(Path2Run), '190603_1039'))=[]; % No vocalization
%  Path2Run(contains(ListOfPaths(Path2Run), '190606_1540'))=[]; % This is not operant but free session, error in choosing the right expe
%  Path2Run(contains(ListOfPaths(Path2Run), '20190712'))=[]; % No neural data
%%
fprintf(1, 'Running result operant bat on %d sessions', length(Path2Run))
for pp=1:length(Path2Run)
    
    Path2ParamFile = ListOfPaths{Path2Run(pp)};
    fprintf(1,'\n\n\n\nRunning result_operant_bat on %d/%d %s\n\n', pp, length(Path2Run),Path2ParamFile)
    
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190202/HoHa_190202_1046_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190131/HoHa_190131_1108_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190130/HoHa_190130_1007_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190129/HoHa_190129_1023_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190120/HoHa_190120_1208_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190119/HoHa_190119_1158_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190118/HoHa_190118_1027_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190117/HoHa_190117_1008_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190116/HoHa_190116_1126_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190122/HoHa_190122_1027_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190123/HoHa_190123_0943_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190124/HoHa_190124_0957_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190125/HoHa_190125_0925_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190128/HoHa_190128_1022_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190203/HoHa_190203_1259_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190204/HoHa_190204_1051_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190205/HoHa_190205_1140_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190206/HoHa_190206_1024_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190207/HoHa_190207_1136_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190208/HoHa_190208_1018_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190211/HoHa_190211_1152_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190212/HoHa_190212_1033_VocTrigger_param.txt';
    % Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190213/HoHa_190213_1101_VocTrigger_param.txt';
%     Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC_HoHa/audio/20190214/HoHa_190214_1130_VocTrigger_param.txt';
%     Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC/LMC_CoEd/audio/20190604/CoEd_190604_1200_VocTrigger_param.txt'; % Needs to point to a reconly param files
%     Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC/LMC_CoEd/audio/20190610/CoEd_190610_0953_VocTrigger_param.txt';
%     Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC/LMC_CoEd/audio/20190607/CoEd_190607_0827_VocTrigger_param.txt';
%     Path2ParamFile = '/Volumes/server_home/users/JulieE/LMC/LMC_CoEd/audio/20190612/CoEd_190612_1030_VocTrigger_param.txt';
%     
     result_operant_bat_what(Path2ParamFile)
    
end

%% RUN audio data and other behavior extraction for the reconly sessions
List2RecOnlyPath = gather_reconly_datapath(BasePath);
%%
Path2RunRecOnly = 1:length(List2RecOnlyPath);
% Path2RunRecOnly = find(contains(List2RecOnlyPath, 'CoEd'));
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '201905'))=[];% No neural data
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '190605_1553'))=[]; % Clock jump, no data to extract
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190703'))=[];% No TTL Pulses
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190709'))=[];% No TTL Pulses
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190619'))=[];% Issue of clock drift for logger 49 and 12

% Path2RunRecOnly = find(contains(List2RecOnlyPath, 'HoHa'));
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190131'))=[];% Allignment issue as of now no neural data extracted
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190202_1400'))=[];% No logger data

%%
for pp= 32:length(Path2RunRecOnly) 

    Path2ParamFile = List2RecOnlyPath{Path2RunRecOnly(pp)};
    fprintf(1,'\n\n\n\nRunning result_reconly_bat on %d/%d %s\n\n',pp, length(Path2RunRecOnly), Path2ParamFile)
    result_reconly_bat(Path2ParamFile)
end

%% Generate the list of paths to gather the data
% BasePath = '/Volumes/server_home/users/JulieE/LMC';
% BasePath = '/Volumes/Julie8T';
%OutputPath = '/Users/elie/Documents/LMCResults';
OutputPath = '/Volumes/JulieE8T/LMCResults';
[ListSSU] = gather_neural_datapath(BasePath);
save(fullfile(OutputPath,'ListSSU.mat'), 'ListSSU')
% Define the path were the data will be saved
% OutputPath = fullfile(BasePath, 'ResultsFiles');
% Files2Run = find(contains(ListSSU,'LMC_CoEd') .* ~contains(ListSSU, '20200109'));
% Files2Run = find(contains(ListSSU,'65701') .* ~contains(ListSSU, '20200109'));
% Files2Run = find(contains(ListSSU,'LMC_HoHa'));
% Files2Run = find(contains(ListSSU,'59882')
%% Sanitary check of neurons
% ( calculate the average spike rate over the...
% whole experiment, measure stability, quality...
% of spike sorting.
fprintf('NEURONS SANITARY CHECK.... ')
Files2Run = 1:length(ListSSU);
% Files2Run = [1:29 87:108];
%  Files2Run=1:488 cells for Co; 123 cells for Ho
for ss=1:length(Files2Run)
    fprintf(1,'File %d/%d: %s\n',ss,length(Files2Run),ListSSU{Files2Run(ss)})
    sanitary_check_perSSfile(ListSSU{Files2Run(ss)}, OutputPath)
end
fprintf(' DONE \n')
% Data for each unit are saved under: sprintf('%s_%s_SS%s_%s-%s.mat', SubjectID, Date,SSQ,TetrodeID,SSID)

%% Sort unit
% Multi-unit SSM should have at least one SNR value above 2
% Single unt SSS should have at least one SNR value above 5 and ISI
% below 0.1%
fprintf('ONLY KEEPING GOOD UNITS.... ')
SSQ_Files2Run = cell(length(Files2Run),1);
for ss=1:length(Files2Run)
    fprintf(1,'File %d/%d\n',ss,length(Files2Run))
    ff=Files2Run(ss);
    [~,FileName] = fileparts(ListSSU{ff});
    Ind_ = strfind(FileName,'_');
    SubjectID = FileName(1:5);
    Date = FileName(7:14);
    TetrodeID = FileName(Ind_(3)-1);
    SSQ = FileName(Ind_(3)+3);
    SSID = FileName((Ind_(4)+1):end);
    Data=load(fullfile(OutputPath,sprintf('%s_%s_SS%s_%s-%s.mat', SubjectID, Date,SSQ,TetrodeID,SSID)));
    if any(Data.QualitySSU.SNR>=5) && (Data.QualitySSU.ISIViolation<=0.1)
        SSQ_Files2Run{ss} = 'SSSU';
    elseif any(Data.QualitySSU.SNR>=2)
        SSQ_Files2Run{ss} = 'SSMU';
    else
        SSQ_Files2Run{ss} = 'NOISE';
    end
end
GoodCellIndices = find(contains(SSQ_Files2Run, 'SS'));
fprintf(' DONE \n')
save('GoodCellIndicesAll.mat','GoodCellIndices','ListSSU','SSQ_Files2Run')

%% Extract the neural data corresponding to the bouts of vocalizations identified
% by voc_localize and voc_localize_operant (run by result_operant_bat.m) for each cell
fprintf(' EXTRACTING NEURAL DATA CORRESPONDING TO VOCALIZATIONS.... \n')
NeuralBuffer = 5000; %duration of the time buffer in
%       ms that should be added before and after the onset and offset time
%       of vocalizations for extracting neural data.

for ss=1:length(GoodCellIndices)
    fprintf(1,'Cell %d/%d: %s \n',ss,length(GoodCellIndices),ListSSU{GoodCellIndices(ss)})
    cut_neuralData_voc_perfile(ListSSU{GoodCellIndices(ss)}, OutputPath,NeuralBuffer)
end
fprintf(' DONE \n')
% Data for each unit and each experimental session are saved as sprintf('%s_%s_%s_SS%s_%s-%s.mat', SubjectID, Date, ExpStartTime,SSQ,TetrodeID,SSID)
% Last run of all dataset on April 29 2022

%% Extract the neural data corresponding to the behaviors identified during the free session
% by get_logger_data_behav (run by result_reconly_bat.m) for each cell
fprintf(' EXTRACTING NEURAL DATA CORRESPONDING TO OTHER BEHAVIORS.... ')


for ss=1:length(GoodCellIndices)
    fprintf(1,'Cell %d/%d\n',ss,length(GoodCellIndices))
    cut_neuralData_behav_perfile(ListSSU{GoodCellIndices(ss)}, OutputPath)
end
fprintf(' DONE \n')
% Data for each unit and each experimental session are saved as or apppended to sprintf('%s_%s_%s_SS%s_%s-%s.mat', SubjectID, Date, ExpStartTime,SSQ, TetrodeID,SSID)


%% Organizing the data as a single file for all behaviors
fprintf(' COMPILING NEURAL DATA .... ')
% turn off warnings for python saving issues
id = 'MATLAB:Python:UnsupportedLoad';
warning('off',id)


for ss=1:length(GoodCellIndices)
    fprintf(1,'Cell %d/%d  %s\n',ss,length(GoodCellIndices), ListSSU{GoodCellIndices(ss)})
    neuralData_compile_perfile(ListSSU{GoodCellIndices(ss)}, OutputPath, NeuralBuffer)
end
warning('on',id)
fprintf(' DONE \n')
% Data for each unit and all experimental session are appended to: sprintf('%s_%s_SS%s_%s-%s.mat', SubjectID, Date,SSQ,TetrodeID,SSID) 
% Last run 04/29/2022 for all dataset

%% Calculating the average spike rate during various types of behaviors including vocalizations
fprintf(' CALCULATING SPIKE RATE CORRESPONDING TO ALL BEHAVIORS.... ')

for ss=1:length(GoodCellIndices)
    fprintf(1,'Cell %d/%d\n',ss,length(GoodCellIndices))
    cal_spikerate_perfile(ListSSU{GoodCellIndices(ss)},OutputPath)
end
fprintf(' DONE \n')
% Data for each unit and all experimental session are appended to: sprintf('%s_%s_SS%s_%s-%s.mat', SubjectID, Date,SSQ,TetrodeID,SSID) 


%% calculate the KDE SPIKE RATE of vocalizations
fprintf(1,' CALCULATING KDE OF THE TIME-VARYING SPIKE RATE CORRESPONDING TO VOCALIZATIONS\n');
Delay = [5000 5000];
for ss=1:length(GoodCellIndices)
    fprintf(1,'Cell %d/%d\n',ss,length(GoodCellIndices))
    cal_kderatevoc_perfile(ListSSU{GoodCellIndices(ss)}, OutputPath,Delay)
end
fprintf(' DONE \n')


%% Plot the average spike rate during various types of behaviors including vocalizations
fprintf(' PLOTING NEURAL DATA (Av RATE) CORRESPONDING TO ALL BEHAVIORS.... ')
UseOldData=1;
TestNames = ["Voc-Free-SelfVsBgd","Self-Free-VocalizingVsQuiet","Voc-Free-SelfVsOthers","Voc-Free-OthersVsBgd","Self-Free-ChewingVsQuiet","Self-Free-LickingVsQuiet"];
SRpValues_all = table('Size',[length(GoodCellIndices) 7], 'VariableTypes',...
    ["logical","double", "double","double","double", "double", "double"],...
    'VariableNames',["Single-Unit" TestNames],...
    'RowNames',ListSSU(GoodCellIndices));
SRcoeffEstimates_all = table('Size',[length(GoodCellIndices) 7], 'VariableTypes',...
    ["logical","double", "double","double","double", "double", "double"],...
    'VariableNames',["Single-Unit" TestNames],...
    'RowNames',ListSSU(GoodCellIndices));

for ss=1:length(GoodCellIndices)
    fprintf(1,'Cell %d/%d\n',ss,length(GoodCellIndices))
    if UseOldData
        [~, DataFile]=fileparts(ListSSU{GoodCellIndices(ss)});
        % Input
        % Get the date of the recording
        Idx_ = strfind(DataFile, '_');
        Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));
        % Get the tetrode ID
        NeuralInputID{1} = DataFile(strfind(DataFile, 'TT')+2);
        % Get the SS ID
        NeuralInputID{2} = DataFile((Idx_(end)+1):end);
        % Get the SS Quality
        NeuralInputID{3} = DataFile(strfind(DataFile, '_SS')+3);
        % Get the subject ID
        SubjectID = DataFile(1:5);
        % Input
        FullDataSetFile = fullfile(OutputPath, sprintf('%s_%s_SS%s_%s-%s.mat', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2}));
        load(FullDataSetFile, 'SpikeRate');
        SRStats = SpikeRate.Stats;

    else
        [SRStats] = plot_av_spikerate_perfile(ListSSU{GoodCellIndices(ss)}, OutputPath);
    end
    SRpValues_all(ss,1) = {strcmp(SSQ_Files2Run{GoodCellIndices(ss)} , 'SSSU')};
    SRcoeffEstimates_all(ss,1) = {strcmp(SSQ_Files2Run{GoodCellIndices(ss)} , 'SSSU')};
    for pp=1:length(TestNames)
        Row = find(strcmp(SRStats.Test, TestNames{pp}));
        if ~isempty(Row)
            SRpValues_all(ss,pp+1) = SRStats(Row,3);
            SRcoeffEstimates_all(ss,pp+1) = SRStats(Row,2);
        else
            SRpValues_all(ss,pp+1) = {nan};
            SRcoeffEstimates_all(ss,pp+1) = {nan};
        end
    end
end
fprintf(' DONE \n')
save(fullfile(OutputPath,'TimeAverageSpikeRate_pValues.mat'), 'SRpValues_all', 'SRcoeffEstimates_all');
%% 

% obtain scatter plot of the coefficient estimates of GLM Poisson for the population data
plim = 0.001;
F1=figure();
tiledlayout(2,1)
nexttile
ScatterMarkerSz = 30;
ColorMU = [0 0.4470 0.7410];
ColorSU = [0.9290 0.6940 0.1250];
BarGraph=nan(4,6);
MU_logical = ~SRpValues_all.("Single-Unit");
NiceNames = {'Vocalizing','Vocalizing vs Quiet','Vocalizing vs Hearing','Hearing', 'Chewing vs Quiet', 'Licking vs Quiet'};
for tt=1:length(TestNames)
    Sig_logical = SRpValues_all.(sprintf(TestNames(tt)))<=plim;
    NSig_logical = SRpValues_all.(sprintf(TestNames(tt)))>plim; % this is to make sure nan values stays 0 indices
    BarGraph(3,tt) = sum(MU_logical.* (~isnan(SRpValues_all.(sprintf(TestNames(tt))))));
    BarGraph(4,tt) = sum((~MU_logical) .* (~isnan(SRpValues_all.(sprintf(TestNames(tt))))));
    BarGraph(1,tt) = sum(MU_logical .* Sig_logical)./BarGraph(3,tt);
    BarGraph(2,tt) = sum((~MU_logical) .* Sig_logical)./BarGraph(4,tt);
    
    swarmchart(tt.*ones(sum(MU_logical.*Sig_logical),1),SRcoeffEstimates_all.(sprintf(TestNames(tt)))(logical(MU_logical.*Sig_logical)),ScatterMarkerSz,ColorMU,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',1)
    hold on
    swarmchart(tt.*ones(sum(MU_logical.*(NSig_logical)),1),SRcoeffEstimates_all.(sprintf(TestNames(tt)))(logical(MU_logical.*NSig_logical)),ScatterMarkerSz,ColorMU,'o','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',1)
    hold on
    swarmchart(tt.*ones(sum(~MU_logical.*(Sig_logical)),1),SRcoeffEstimates_all.(sprintf(TestNames(tt)))(logical(~MU_logical.*Sig_logical)),ScatterMarkerSz,ColorSU,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',1)
    hold on
    swarmchart(tt.*ones(sum(~MU_logical.*(NSig_logical)),1),SRcoeffEstimates_all.(sprintf(TestNames(tt)))(logical(~MU_logical.*NSig_logical)),ScatterMarkerSz,ColorSU,'o','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',1)
    text(tt-.25,10^-1.3,sprintf('%d MU',BarGraph(3,tt)), 'Color',ColorMU)
    text(tt-.25,10^-1.5,sprintf('%d SU',BarGraph(4,tt)), 'Color',ColorSU)
    if tt==1
        legend(sprintf('MU p<=%.3f',plim),sprintf('MU p>%.3f',plim), sprintf('SU p<=%.3f',plim),sprintf('SU p>%.3f',plim),  'boxoff', 'Location','north', 'Orientation', 'horizontal')
        legend('AutoUpdate','off')
    end
end
legend('boxoff')
F1.Children.Children(2).XTick = 1:length(TestNames);
F1.Children.Children(2).XTickLabel = NiceNames;
F1.Children.Children(2).XTickLabelRotation = 25;
YLIM = F1.Children.Children(2).YLim;
F1.Children.Children(2).YLim = [10^-1.7 10^2.2];
F1.Children.Children(2).YScale = 'log';
ylabel('Coefficient estimates for each unit')
% title('Poisson GLM on time average spike rate')


nexttile
X = categorical(NiceNames);
X = reordercats(X,NiceNames);
B=bar(X,BarGraph(1:2,:),'EdgeColor','k')
ylabel(sprintf('Proportion of significant units at p<=%.3f',plim))

B(1).Parent.YLim = [0 1];
B(1).FaceColor = ColorMU;
B(2).FaceColor = ColorSU;
for bb=1:size(BarGraph,2)
    text(bb-.25,0.97,sprintf('%d MU',BarGraph(3,bb)), 'Color', ColorMU)
    text(bb-.25,0.93,sprintf('%d SU',BarGraph(4,bb)), 'Color', ColorSU)
end
% title('Poisson GLM on time average spike rate')
suplabel('Poisson GLM on time average spike rate', 't')

% find the single unit with max estimate during vocalization
Sig_logical = SRpValues_all.(sprintf(TestNames(1)))<=plim;
[~,IndSUdesc] = sort(SRcoeffEstimates_all.(sprintf(TestNames(1)))(logical(~MU_logical.*Sig_logical)), 'descend');
IndSUMax = find(SRcoeffEstimates_all.(sprintf(TestNames(1)))(logical(~MU_logical.*Sig_logical)) == max(SRcoeffEstimates_all.(sprintf(TestNames(1)))(logical(~MU_logical.*Sig_logical))));
ListSSU_SUSig = ListSSU(GoodCellIndices(logical(~MU_logical.*Sig_logical)));
fprintf(1,'the ID of the SU with max estimate for vocalizing vs Background is %s', ListSSU_SUSig{IndSUMax})

% The plot is saved under OutputPath as sprintf('%s_%s_%s_SS%s_%s-%s_MeanRateScatter.pdf', SubjectID, SSQ,TetrodeID,SSID))
%% Plot rasters for vocalizations
fprintf(1,' RASTER PLOTS (AND KDE) of NEURAL DATA CORRESPONDING TO VOCALIZATIONS\n');
Delay = [5000 5000];
PlotDyn = 0; %Set to 1 to plot dnamic plots
DurOrd = 0; % set to 1 to order neural responses by increasing vocalization duration
for ss=1:length(GoodCellIndices)
    fprintf(1,'Cell %d/%d\n',ss,length(GoodCellIndices))
    plot_rastervoc_perfile(ListSSU{Files2Run(GoodCellIndices(ss))}, OutputPath, Delay, PlotDyn, DurOrd)
    close all
end
fprintf(' DONE \n')

%% Plot the time varying rate of vocalizations
fprintf(1,' PLOTING KDE OF THE TIME-VARYING SPIKE RATE CORRESPONDING TO VOCALIZATIONS\n');

for ss=475:length(GoodCellIndices)
    fprintf(1,'Cell %d/%d\n',ss,length(GoodCellIndices))
    plot_kderatevoc_perfile(ListSSU{Files2Run(GoodCellIndices(ss))}, OutputPath, Delay)
end
fprintf(' DONE \n')

    









%% Get the path to audio data for operant conditioning experiment
    [AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
    Date = DataFile(6:11);
    fprintf(1,'Working on %s\n', Date)
    ExpStartTime = DataFile(13:16);
    % Set the path to logger data
    Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'logger',['20' Date]);
    
    % Set the time buffer before vocalizations onset
    BufferBeforeOnset = 200; %ms
    fprintf(' EXTRACTING NEURAL DATA CORRESPONDING TO VOCALIZATIONS \n')
    FlagsExtr = [0 0 1 1]; % FlagsExtr(1)= Raw data, FlagsExtr(2) = LFP, FlagsExtr(3) = Tetrodes, FlagsExtr(4) = single units
    DenoiseT=1;
    Rthreshold = [0.92 0.94 0.96 0.98];
    cut_neuralData_voc(Logger_dir,Date, ExpStartTime,FlagsExtr,BufferBeforeOnset,DenoiseT,Rthreshold);
    
    %% Plot PSTH of the bats hearing or producing a vocalization during the operant conditioning
    % Find the ID of the Neural loggers and corresponding audiologger for each implanted bat
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:P200','basic');
    RowData = find((cell2mat(RecTableData(2:end,1))== str2double(Date))) +1;
    DataInfo = RecTableData(RowData,:);
    Header = RecTableData(1,:);
    BatIDCol = find(contains(Header, 'Bat'));
    NLCol = find(contains(Header, 'NL'));
    ALCol = find(contains(Header, 'AL-throat'));
    NL_ID = cell2mat(DataInfo(NLCol));
    for nl=1:length(NL_ID)
        NeuroLoggerID = ['Logger' num2str(NL_ID(nl))];
        AL_ID = DataInfo{ALCol(find(ALCol<NLCol(nl),1,'last'))};
        AudioLoggerID = ['Logger' num2str(AL_ID)];
        Flags=[0 1];% Flags = whether to calculate PSTH of Tetrode (Flags(1)=1) and/or Single units
        % (Flags(2)=1))
        KDE_Cal = 1;
        PLOT=1; % set to 1 to plot the results, 0 to just return data
        fprintf(1,' PSTH of NEURAL DATA CORRESPONDING TO VOCALIZATIONS DURING OPERANT %s\n', NeuroLoggerID)
%         for rr=1:length(Rthreshold)
%             [SpikeTimesVoc.(NeuroLoggerID)]=plot_psth_voc(Logger_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, BufferBeforeOnset, KDE_Cal,PLOT,rr);
%         end
        [SpikeTimesVoc.(NeuroLoggerID)]=plot_psth_voc(Logger_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, BufferBeforeOnset, KDE_Cal,PLOT);
        [SpikeTimesVocOff.(NeuroLoggerID)]=plot_psth_voc_off(Logger_dir, Date, ExpStartTime, AudioLoggerID, NeuroLoggerID, Flags, BufferBeforeOnset, KDE_Cal,PLOT);
    end
    save(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, BufferBeforeOnset)), 'SpikeTimesVoc','-append')
    save(fullfile(Logger_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, BufferBeforeOnset)), 'SpikeTimesVocOff','-append')
    close all
    %     pause()
    
     %% Extract data of the bats doing other actions during the free behavior session (RecOnly)
    fprintf(' EXTRACTING ONSET/OFFSET TIMES OF OTHER BEHAVIORS DURING FREE SESSION \n')
    RecOnlySession = dir(fullfile(AudioDataPath, '*RecOnly_events.txt'));
    if isempty(RecOnlySession)
        fprintf(1,'No free interaction session on that day!\n')
    else
        if length(RecOnlySession)>1
            fprintf(1, 'Several RecOnly session were done on that day:\n')
            for ss=1:length(RecOnlySession)
                fprintf(1, '%d. %s\n', ss, RecOnlySession(ss).name);
            end
            Inputss = input('Your choice:\n');
            RecOnlySession = RecOnlySession(Inputss);
        end
        Date = RecOnlySession.name(6:11);
        ExpStartTime = RecOnlySession.name(13:16);
        % extract the time onset/offset of behaviors
        get_logger_data_behav(AudioDataPath, Logger_dir, Date, ExpStartTime)
        %
        %% Extract the neural data corresponding to other actions during the free behavior session (RecOnly)
        fprintf(' EXTRACTING NEURAL DATA CORRESPONDING TO OTHER BEHAVIORS DURING FREE SESSION \n')
        FlagsExtr = [0 0 1 1]; % FlagsExtr(1)= Raw data, FlagsExtr(2) = LFP, FlagsExtr(3) = Tetrodes, FlagsExtr(4) = single units
        BufferBeforeBehavOnset = 0;
        RecOnlySession = dir(fullfile(AudioDataPath, '*RecOnly_events.txt'));
        if length(RecOnlySession)>1
            fprintf(1, 'Several RecOnly session were done on that day:\n')
            for ss=1:length(RecOnlySession)
                fprintf(1, '%d. %s\n', ss, RecOnlySession(ss).name);
            end
            Inputss = input('Your choice:\n');
            RecOnlySession = RecOnlySession(Inputss);
        end
        Date = RecOnlySession.name(6:11);
        ExpStartTime = RecOnlySession.name(13:16);
        cut_neuralData_behav(Logger_dir,Date, ExpStartTime,FlagsExtr,BufferBeforeBehavOnset);
        %
        %% Plot PSTH of the bats doing other actions during free socialization!
        RecOnlySession = dir(fullfile(AudioDataPath, '*RecOnly_events.txt'));
        if length(RecOnlySession)>1
            fprintf(1, 'Several RecOnly session were done on that day:\n')
            for ss=1:length(RecOnlySession)
                fprintf(1, '%d. %s\n', ss, RecOnlySession(ss).name);
            end
            Inputss = input('Your choice:\n');
            RecOnlySession = RecOnlySession(Inputss);
        end
        Date = RecOnlySession.name(6:11);
        ExpStartTime = RecOnlySession.name(13:16);
        % Find the ID of the Neural loggers and corresponding ID tag of each implanted bat
        MaxDur = 700; % duration by which each long behavioral sequence should be cut often set at 600 but for 190118, 700 is better
        [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:P200','basic');
        RowData = find((cell2mat(RecTableData(2:end,1))== str2double(Date))) +1;
        DataInfo = RecTableData(RowData,:);
        Header = RecTableData(1,:);
        BatIDCol = find(contains(Header, 'Bat'));
        NLCol = find(contains(Header, 'NL'));
        NL_ID = cell2mat(DataInfo(NLCol));
        for nl=1:length(NL_ID)
            NeuroLoggerID = ['Logger' num2str(NL_ID(nl))];
            Bat_ID = DataInfo{BatIDCol(find(BatIDCol<NLCol(nl),1,'last'))};
            Flags=[1 1];% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
            % (Flags(2)=1))
            PLOT = 0; % set to 1 to plot the results, 0 to just return data
            KDE_Cal = 1;
            fprintf(1,' PSTH of NEURAL DATA CORRESPONDING TO OTHER BEHAVIORS DURING FREE SOCIALIZATION %s \n', NeuroLoggerID)
            [SpikeTimesBehav.(NeuroLoggerID)]= plot_psth_behav(Logger_dir, Date, ExpStartTime, NeuroLoggerID,Bat_ID, Flags, MaxDur, KDE_Cal, PLOT);
        end
        save(fullfile(Logger_dir, sprintf('%s_%s_BehavExtractData_%d.mat', Date, ExpStartTime,BufferBeforeBehavOnset)),'SpikeTimesBehav','-append');
        close all
    end
%     
    %% Plot one PSTH per unit with all actions
    if isempty(RecOnlySession)
            fprintf(1,'No free interaction session on that day!\n')
            for nl=1:length(NL_ID)
                NeuroLoggerID = ['Logger' num2str(NL_ID(nl))];
                SpikeTimesBehav.(NeuroLoggerID) = [];
            end
    else
        load(fullfile(Logger_dir, sprintf('%s_%s_BehavExtractData_%d.mat', Date, ExpStartTime,BufferBeforeBehavOnset)) , 'SpikeTimesBehav');
    end
    for nl=1:length(NL_ID)
        NeuroLoggerID = ['Logger' num2str(NL_ID(nl))];
        Bat_ID = DataInfo{BatIDCol(find(BatIDCol<NLCol(nl),1,'last'))};
        Flags=[1 1];% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
        % (Flags(2)=1))
        KDE_Cal = 1;
        
        fprintf(1,' PSTH of NEURAL DATA CORRESPONDING TO BEHAVIORS DURING FREE SOCIALIZATION AND VOCAL ACTIVITY DURING OPERANT CONDITIONING %s\n', NeuroLoggerID)
        plot_psth_voc_and_behav(SpikeTimesBehav.(NeuroLoggerID),SpikeTimesVoc.(NeuroLoggerID),Logger_dir,Date, NeuroLoggerID,Flags, BufferBeforeOnset,MaxDur, KDE_Cal);
        plot_psth_voc_and_behav_off(SpikeTimesBehav.(NeuroLoggerID),SpikeTimesVocOff.(NeuroLoggerID),Logger_dir,Date, NeuroLoggerID,Flags, BufferBeforeOnset,MaxDur, KDE_Cal);
    end

    
%%   INTERNAL FUNCTIONS 
function [ListSSU] = gather_neural_datapath(BasePath)
fprintf(1,'*** Gathering paths to spike sorted units ***')
ListSSU = cell(10^3,1); % initialize the list to 1000
ExpFolders = dir(fullfile(BasePath,'LMC*'));
NSSU = 0; % counter for single units
for ee=1:length(ExpFolders)
    fprintf(1, '\n  -> Looking into  %s...\n ', fullfile(ExpFolders(ee).folder,ExpFolders(ee).name))
    DateFolders = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name, 'logger','20*'));
    for dd=1:length(DateFolders)
        fprintf(1, '   %s\n', DateFolders(dd).name);
        LoggerFolders = dir(fullfile(DateFolders(dd).folder, DateFolders(dd).name,'Logger*'));
        for ll = 1:length(LoggerFolders)
            SSFiles = dir(fullfile(LoggerFolders(ll).folder, LoggerFolders(ll).name, 'extracted_data', '*_TT*_SSS*.mat'));
            if ~isempty(SSFiles)
                for ssf=1:length(SSFiles)
                    NSSU = NSSU +1;
                    ListSSU{NSSU} = fullfile(SSFiles(ssf).folder, SSFiles(ssf).name);
                end
            end
            SSFiles = dir(fullfile(LoggerFolders(ll).folder, LoggerFolders(ll).name, 'extracted_data', '*_TT*_SSM*.mat'));
            if ~isempty(SSFiles)
                for ssf=1:length(SSFiles)
                    NSSU = NSSU +1;
                    ListSSU{NSSU} = fullfile(SSFiles(ssf).folder, SSFiles(ssf).name);
                end
            end
            SSFiles = dir(fullfile(LoggerFolders(ll).folder, LoggerFolders(ll).name, 'extracted_data', '*_TT*_SSU*.mat'));
            if ~isempty(SSFiles)
                for ssf=1:length(SSFiles)
                    NSSU = NSSU +1;
                    ListSSU{NSSU} = fullfile(SSFiles(ssf).folder, SSFiles(ssf).name);
                end
            end
        end
    end
end
ListSSU = ListSSU(1:NSSU);
fprintf(1, '\n Files from %d single units have been retrieved\n', NSSU);
end


function [List2ParamPath] = gather_operant_datapath(BasePath)
fprintf(1,'*** Gathering paths to audio operant data ***')
List2ParamPath = cell(10^3,1); % initialize the list to 1000
ExpFolders = dir(fullfile(BasePath,'LMC*'));
NF = 0; % counter for single files
for ee=1:length(ExpFolders)
    fprintf(1, '\n  -> Looking into  %s...\n ', fullfile(ExpFolders(ee).folder,ExpFolders(ee).name))
    DateFolders = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name, 'audio','20*'));
    for dd=1:length(DateFolders)
        fprintf(1, '   %s\n', DateFolders(dd).name);
        AudioParamFiles = dir(fullfile(DateFolders(dd).folder, DateFolders(dd).name,'*VocTrigger_param.txt'));
        for ll = 1:length(AudioParamFiles)
            NF = NF +1;
            List2ParamPath{NF} = fullfile(AudioParamFiles(ll).folder, AudioParamFiles(ll).name);
        end
    end
end
List2ParamPath = List2ParamPath(1:NF);
fprintf(1, '\n Files from %d sessions or operant conditioning have been retrieved\n', NF);
end

function [List2ParamPath] = gather_reconly_datapath(BasePath)
fprintf(1,'*** Gathering paths to audio reconly data ***')
List2ParamPath = cell(10^3,1); % initialize the list to 1000
ExpFolders = dir(fullfile(BasePath,'LMC*'));
NF = 0; % counter for single files
for ee=1:length(ExpFolders)
    fprintf(1, '\n  -> Looking into  %s...\n ', fullfile(ExpFolders(ee).folder,ExpFolders(ee).name))
    DateFolders = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name, 'audio','20*'));
    for dd=1:length(DateFolders)
        fprintf(1, '   %s\n', DateFolders(dd).name);
        AudioParamFiles = dir(fullfile(DateFolders(dd).folder, DateFolders(dd).name,'*RecOnly_param.txt'));
        for ll = 1:length(AudioParamFiles)
            NF = NF +1;
            List2ParamPath{NF} = fullfile(AudioParamFiles(ll).folder, AudioParamFiles(ll).name);
        end
    end
end
List2ParamPath = List2ParamPath(1:NF);
fprintf(1, '\n Files from %d sessions or free recording sessions have been retrieved\n', NF);
end