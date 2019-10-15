function cal_kderatevoc_perfile(InputDataFile)
%%
if nargin<3
    Delay=[3000 200];% in ms
end

Bin_ms = 1; % size of the KDE binning
%Response_samprate = 1/Bin_ms;% Sampling rate of the KDE in kHz

[~, DataFile]=fileparts(InputDataFile);
% Input
% Get the date of the recording
Idx_ = strfind(DataFile, '_');
Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));

% Get the tetrode ID
NeuralInputID{1} = DataFile(strfind(DataFile, 'TT')+2);
% Get the SS ID
NeuralInputID{2} = DataFile((Idx_(end)+1):end);

% Get the subject ID
SubjectID = DataFile(1:5);

% Input
FileNameBase = sprintf('%s_%s_SSU%s-%s', SubjectID, Date,NeuralInputID{1},NeuralInputID{2});
FullDataSetFile = fullfile(OutputPath, sprintf('%s.mat', FileNameBase));
Data=load(FullDataSetFile);

% Indices of vocalization renditions produced by subject ('self') or not
IndVocP = find(contains(Data.What, 'Voc').*contains(Data.Who, 'self'));
IndVocH = find(contains(Data.What, 'Voc').*~contains(Data.Who, 'self'));
IndVocD = find(((Data.DelayBefore>=Delay(1)) + contains(Data.VocRank, 'first')).*((Data.DelayAfter>=Delay(2))+ contains(Data.VocRank, 'end')));
IndVocO = find(contains(Data.ExpType, 'O'));
IndVocF = find(contains(Data.ExpType, 'F'));
IndVocPD = intersect(IndVocP, IndVocD);
IndVocPDO = intersect(IndVocPD, IndVocO);
IndVocPDF = intersect(IndVocPD, IndVocF);
IndVocHD = intersect(IndVocH, IndVocD);
IndVocHDO = intersect(IndVocHD, IndVocO);
IndVocHDF = intersect(IndVocHD, IndVocF);
IndTr = find(contains(Data.What, 'VocTr'));
IndBa = find(contains(Data.What, 'VocBa'));

%% KDE time-varying rate alligned to vocalization production onset/offset
if ~isempty(IndVocPD) && ~isempty(IndVocPDO) && ~isempty(IndVocPDF)
    [KDE.SelfVocAll] = kderate(Data.SpikesArrivalTimes_Behav(IndVocPD),Data.Duration(IndVocPD),Delay,Bin_ms);
    IndTrPD = intersect(IndTr, IndVocPD);
    if length(IndTrPD)>MinNumCall
        [KDE.SelfVocTrAll] = kderate(Data.SpikesArrivalTimes_Behav(IndTrPD),Data.Duration(IndTrPD),Delay,Bin_ms);
    end
    IndBaPD = intersect(IndBa, IndVocPD);
    if length(IndBaPD)>MicNumCall
        [KDE.SelfVocBaAll] = kderate(Data.SpikesArrivalTimes_Behav(IndBaPD),Data.Duration(IndBaPD),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization perception onset/offset
if ~isempty(IndVocHD) && ~isempty(IndVocHDO) && ~isempty(IndVocHDF)
    [KDE.OthersVocAll] = kderate(Data.SpikesArrivalTimes_Behav(IndVocHD),Data.Duration(IndVocHD),Delay,Bin_ms);
    IndTrHD = intersect(IndTr, IndVocHD);
    if length(IndTrHD)>MinNumCall
        [KDE.OthersVocTrAll] = kderate(Data.SpikesArrivalTimes_Behav(IndTrHD),Data.Duration(IndTrHD),Delay,Bin_ms);
    end
    IndBaHD = intersect(IndBa, IndVocHD);
    if length(IndBaHD)>MinNumCall
        [KDE.OthersVocBaAll] = kderate(Data.SpikesArrivalTimes_Behav(IndBaHD),Data.Duration(IndBaHD),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization production onset/offset during Operant conditioning
if ~isempty(IndVocPDO)
    [KDE.SelfVocOp] = kderate(Data.SpikesArrivalTimes_Behav(IndVocPDO),Data.Duration(IndVocPDO),Delay,Bin_ms);
    IndTrPDO = intersect(IndTr, IndVocPDO);
    if length(IndTrPDO)>MinNumCall
        [KDE.SelfVocTrOp] = kderate(Data.SpikesArrivalTimes_Behav(IndTrPDO),Data.Duration(IndTrPDO),Delay,Bin_ms);
    end
    IndBaPDO = intersect(IndBa, IndVocPDO);
    if length(IndBaPDO)>MinNumCall
        [KDE.SelfVocBaOp] = kderate(Data.SpikesArrivalTimes_Behav(IndBaPDO),Data.Duration(IndBaPDO),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization perception onset/offset during operant conditioning
if ~isempty(IndVocHDO)
    [KDE.OthersVocOp] = kderate(Data.SpikesArrivalTimes_Behav(IndVocHDO),Data.Duration(IndVocHDO),Delay,Bin_ms);
    IndTrHDO = intersect(IndTr, IndVocHDO);
    if length(IndTrHDO)>MinNumCall
        [KDE.OthersVocTrOp] = kderate(Data.SpikesArrivalTimes_Behav(IndTrHDO),Data.Duration(IndTrHDO),Delay,Bin_ms);
    end
    IndBaHDO = intersect(IndBa, IndVocHDO);
    if length(IndBaPDO)>MinNumCall
        [KDE.OthersVocBaOp] = kderate(Data.SpikesArrivalTimes_Behav(IndBaHDO),Data.Duration(IndBaHDO),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization production onset/offset during Free session
if ~isempty(IndVocPDF)
    [KDE.SelfVocFr] = kderate(Data.SpikesArrivalTimes_Behav(IndVocPDF),Data.Duration(IndVocPDF),Delay,Bin_ms);
    IndTrPDF = intersect(IndTr, IndVocPDF);
    if length(IndTrPDF)>MinNumCall
        [KDE.SelfVocTrFr] = kderate(Data.SpikesArrivalTimes_Behav(IndTrPDF),Data.Duration(IndTrPDF),Delay,Bin_ms);
    end
    IndBaPDF = intersect(IndBa, IndVocPDF);
    if length(IndBaPDF)>MinNumCall
        [KDE.SelfVocBaFr] = kderate(Data.SpikesArrivalTimes_Behav(IndBaPDF),Data.Duration(IndBaPDF),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization perception onset/offset during Free session
if ~isempty(IndVocHDF)
    [KDE.OthersVocFr] = kderate(Data.SpikesArrivalTimes_Behav(IndVocHDF),Data.Duration(IndVocHDF),Delay,Bin_ms);
    IndTrHDF = intersect(IndTr, IndVocHDF);
    if length(IndTrHDF)>MinNumCall
        [KDE.OthersVocTrFr] = kderate(Data.SpikesArrivalTimes_Behav(IndTrHDF),Data.Duration(IndTrHDF),Delay,Bin_ms);
    end
    IndBaHDF = intersect(IndBa, IndVocHDF);
    if length(IndBaPDO)>MinNumCall
        [KDE.OthersVocBaFr] = kderate(Data.SpikesArrivalTimes_Behav(IndBaHDF),Data.Duration(IndBaHDF),Delay,Bin_ms);
    end
end



% save data
save(FullDataSetFile, 'KDE', '-append')


%% INTERNAL FUNCTION
    function [OUT] = kderate(SpikesArrivalTimes,Duration,Delay, Bin_ms)
        Response_samprate = 1/Bin_ms;% Sampling rate of the KDE in kHz
        % Find the number of events for each time window
        t=-Delay(1): Bin_ms : round((max(Duration) + Delay(2))/Bin_ms)*Bin_ms;
        Weight = zeros(1,length(t));
        for vv=1:length(Duration)
            t_local = -Delay(1): Bin_ms : (Duration(vv) + Delay(2));
            Weight(1:length(t_local))=Weight(1:length(t_local))+1;
        end
        % calculated KDE
        AllSpikes_local = cell2mat(SpikesArrivalTimes);
        [Kde,T,Error] = kde_wrapper(AllSpikes_local,t,Response_samprate,sum(~isnan(PSTH_local)));
        OUT = [Kde;T;Error];
    end
end
