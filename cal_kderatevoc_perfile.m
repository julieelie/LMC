function cal_kderatevoc_perfile(InputDataFile,OutputPath, Delay)
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
% Get the SS quality
NeuralInputID{3} = DataFile(strfind(DataFile, '_SS')+3);

% Get the subject ID
SubjectID = DataFile(1:5);

% Input
FileNameBase = sprintf('%s_%s_SS%s_%s-%s', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2});
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
    [KDE_onset.SelfVocAll] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndVocPD),Data.Duration(IndVocPD),Delay,Bin_ms);
    [KDE_offset.SelfVocAll] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndVocPD),Data.Duration(IndVocPD),Delay,Bin_ms);
    
    IndTrPD = intersect(IndTr, IndVocPD);
    if length(IndTrPD)>MinNumCall
        [KDE_onset.SelfTrAll] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndTrPD),Data.Duration(IndTrPD),Delay,Bin_ms);
        [KDE_offset.SelfTrAll] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndTrPD),Data.Duration(IndTrPD),Delay,Bin_ms);
    
    end
    IndBaPD = intersect(IndBa, IndVocPD);
    if length(IndBaPD)>MicNumCall
        [KDE_onset.SelfBaAll] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndBaPD),Data.Duration(IndBaPD),Delay,Bin_ms);
        [KDE_offset.SelfBaAll] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndBaPD),Data.Duration(IndBaPD),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization perception onset/offset
if ~isempty(IndVocHD) && ~isempty(IndVocHDO) && ~isempty(IndVocHDF)
    [KDE_onset.OthersVocAll] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndVocHD),Data.Duration(IndVocHD),Delay,Bin_ms);
    [KDE_offset.OthersVocAll] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndVocHD),Data.Duration(IndVocHD),Delay,Bin_ms);
    IndTrHD = intersect(IndTr, IndVocHD);
    if length(IndTrHD)>MinNumCall
        [KDE_onset.OthersTrAll] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndTrHD),Data.Duration(IndTrHD),Delay,Bin_ms);
        [KDE_offset.OthersTrAll] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndTrHD),Data.Duration(IndTrHD),Delay,Bin_ms);
    end
    IndBaHD = intersect(IndBa, IndVocHD);
    if length(IndBaHD)>MinNumCall
        [KDE_onset.OthersBaAll] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndBaHD),Data.Duration(IndBaHD),Delay,Bin_ms);
        [KDE_offset.OthersBaAll] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndBaHD),Data.Duration(IndBaHD),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization production onset/offset during Operant conditioning
if ~isempty(IndVocPDO)
    [KDE_onset.SelfVocOp] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndVocPDO),Data.Duration(IndVocPDO),Delay,Bin_ms);
    [KDE_offset.SelfVocOp] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndVocPDO),Data.Duration(IndVocPDO),Delay,Bin_ms);
    IndTrPDO = intersect(IndTr, IndVocPDO);
    if length(IndTrPDO)>MinNumCall
        [KDE_onset.SelfTrOp] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndTrPDO),Data.Duration(IndTrPDO),Delay,Bin_ms);
        [KDE_offset.SelfTrOp] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndTrPDO),Data.Duration(IndTrPDO),Delay,Bin_ms);
    end
    IndBaPDO = intersect(IndBa, IndVocPDO);
    if length(IndBaPDO)>MinNumCall
        [KDE_onset.SelfBaOp] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndBaPDO),Data.Duration(IndBaPDO),Delay,Bin_ms);
        [KDE_offset.SelfBaOp] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndBaPDO),Data.Duration(IndBaPDO),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization perception onset/offset during operant conditioning
if ~isempty(IndVocHDO)
    [KDE_onset.OthersVocOp] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndVocHDO),Data.Duration(IndVocHDO),Delay,Bin_ms);
    [KDE_offset.OthersVocOp] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndVocHDO),Data.Duration(IndVocHDO),Delay,Bin_ms);
    IndTrHDO = intersect(IndTr, IndVocHDO);
    if length(IndTrHDO)>MinNumCall
        [KDE_onset.OthersTrOp] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndTrHDO),Data.Duration(IndTrHDO),Delay,Bin_ms);
        [KDE_offset.OthersTrOp] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndTrHDO),Data.Duration(IndTrHDO),Delay,Bin_ms);
    end
    IndBaHDO = intersect(IndBa, IndVocHDO);
    if length(IndBaPDO)>MinNumCall
        [KDE_onset.OthersBaOp] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndBaHDO),Data.Duration(IndBaHDO),Delay,Bin_ms);
        [KDE_offset.OthersBaOp] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndBaHDO),Data.Duration(IndBaHDO),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization production onset/offset during Free session
if ~isempty(IndVocPDF)
    [KDE_onset.SelfVocFr] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndVocPDF),Data.Duration(IndVocPDF),Delay,Bin_ms);
    [KDE_offset.SelfVocFr] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndVocPDF),Data.Duration(IndVocPDF),Delay,Bin_ms);
    IndTrPDF = intersect(IndTr, IndVocPDF);
    if length(IndTrPDF)>MinNumCall
        [KDE_onset.SelfTrFr] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndTrPDF),Data.Duration(IndTrPDF),Delay,Bin_ms);
        [KDE_offset.SelfTrFr] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndTrPDF),Data.Duration(IndTrPDF),Delay,Bin_ms);
    end
    IndBaPDF = intersect(IndBa, IndVocPDF);
    if length(IndBaPDF)>MinNumCall
        [KDE_onset.SelfBaFr] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndBaPDF),Data.Duration(IndBaPDF),Delay,Bin_ms);
        [KDE_offset.SelfBaFr] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndBaPDF),Data.Duration(IndBaPDF),Delay,Bin_ms);
    end
end

%% KDE time-varying rate alligned to vocalization perception onset/offset during Free session
if ~isempty(IndVocHDF)
    [KDE_onset.OthersVocFr] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndVocHDF),Data.Duration(IndVocHDF),Delay,Bin_ms);
    [KDE_offset.OthersVocFr] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndVocHDF),Data.Duration(IndVocHDF),Delay,Bin_ms);
    IndTrHDF = intersect(IndTr, IndVocHDF);
    if length(IndTrHDF)>MinNumCall
        [KDE_onset.OthersTrFr] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndTrHDF),Data.Duration(IndTrHDF),Delay,Bin_ms);
        [KDE_offset.OthersTrFr] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndTrHDF),Data.Duration(IndTrHDF),Delay,Bin_ms);
    end
    IndBaHDF = intersect(IndBa, IndVocHDF);
    if length(IndBaPDO)>MinNumCall
        [KDE_onset.OthersBaFr] = kderate_onset(Data.SpikesArrivalTimes_Behav(IndBaHDF),Data.Duration(IndBaHDF),Delay,Bin_ms);
        [KDE_offset.OthersBaFr] = kderate_offset(Data.SpikesArrivalTimes_Behav(IndBaHDF),Data.Duration(IndBaHDF),Delay,Bin_ms);
    end
end



% save data
save(FullDataSetFile, 'KDE_onset','KDE_offset', '-append')


%% INTERNAL FUNCTION
    function [OUT] = kderate_onset(SpikesArrivalTimes,Duration,Delay, Bin_ms)
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
        [Kde,T,Error] = kde_wrapper(AllSpikes_local,t,Response_samprate,Weight); % Calculate the kde in spike /ms
        OUT = [Kde*10^3;T;Error*10^3];% save the kde in Hz (spike /s) and timeline in ms
    end

    function [OUT] = kderate_offset(SpikesArrivalTimes,Duration,Delay, Bin_ms)
        Response_samprate = 1/Bin_ms;% Sampling rate of the KDE in kHz
        % Find the number of events for each time window and center spikes
        % to vocalization offset
        t = -round((max(Duration) + Delay(2))/Bin_ms)*Bin_ms : Bin_ms : Delay(1);
        Weight = zeros(1,length(t));
        for vv=1:length(Duration)
            t_local = -(Duration(vv) + Delay(2)) : Bin_ms : Delay(1);
            Weight((length(t)-length(t_local)+1):length(t))=Weight((length(t)-length(t_local)+1):length(t))+1;
            SpikesArrivalTimes{vv} = SpikesArrivalTimes{vv}-Duration(vv);
        end
        % calculated KDE
        AllSpikes_local = cell2mat(SpikesArrivalTimes);
        [Kde,T,Error] = kde_wrapper(AllSpikes_local,t,Response_samprate,Weight); % Calculate the kde in spike /ms
        OUT = [Kde*10^3;T;Error*10^3];% save the kde in Hz (spike /s) and timeline in ms
    end
end
