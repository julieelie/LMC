function cal_spikerate_perfile(InputDataFile,OutputPath)
%%
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

% Output
FullDataSetFile = fullfile(OutputPath, sprintf('%s_%s_SSU%s-%s.mat', SubjectID, Date,NeuralInputID{1},NeuralInputID{2}));
Data=load(FullDataSetFile);

%% Data pointers
% Indices of vocalization renditions produced by subject ('self') or not
IndVocS = find(contains(Data.What, 'Voc').*contains(Data.Who, 'self'));
IndVocO = find(contains(Data.What, 'Voc').*~contains(Data.Who, 'self'));
IndVoc = find(contains(Data.What, 'Voc'));

% Indices of non-vocal behavior
IndNVBehav = ~contains(Data.What, 'Voc');
% Find the list of other behaviors extracted
BehavType = unique(Data.What(IndNVBehav));

%% Calculate mean rate
% Calculate the spike rate in Hz of self calls
if ~isempty(IndVocS)
    SelfCall_rate=nan(length(IndVocS),2);
    SelfCall_exptype=cell(length(IndVocS),1);
    for oo=1:length(IndVocS)
        SAT = Data.SpikesArrivalTimes_Behav{IndVocS(oo)};
        NumSpikes = sum((SAT>=0).*(SAT<(Data.Duration(IndVocS(oo)))));
        SelfCall_rate(oo,1) = NumSpikes./Data.Duration(IndVocS(oo)).*10^3;
        SAT = Data.SpikesArrivalTimes_Baseline{IndVocS(oo)};
        NumSpikes = sum((SAT>=0).*(SAT<(Data.BSLDuration(IndVocS(oo)))));
        SelfCall_rate(oo,2) = NumSpikes./Data.BSLDuration(IndVocS(oo)).*10^3;
        SelfCall_exptype{oo} = Data.ExpType{IndVocS(oo)};
    end
else
    SelfCall_rate = nan(1,2);
    SelfCall_exptype = {};
end

% Calculate the spike rate in Hz of others' calls
if ~isempty(IndVocO)
    OthersCall_rate=nan(length(IndVocO),2);
    OthersCall_exptype=cell(length(IndVocO),1);
    for oo=1:length(IndVocO)
        SAT = Data.SpikesArrivalTimes_Behav{IndVocO(oo)};
        NumSpikes = sum((SAT>=0).*(SAT<(Data.Duration(IndVocO(oo)))));
        OthersCall_rate(oo,1) = NumSpikes./Data.Duration(IndVocO(oo)).*10^3;
        SAT = Data.SpikesArrivalTimes_Baseline{IndVocO(oo)};
        NumSpikes = sum((SAT>=0).*(SAT<(Data.BSLDuration(IndVocO(oo)))));
        OthersCall_rate(oo,2) = NumSpikes./Data.BSLDuration(IndVocO(oo)).*10^3;
        OthersCall_exptype{oo} = Data.ExpType{IndVocO(oo)};
    end
else
    OthersCall_rate = nan(1,2);
    OthersCall_exptype = {};
end

% Calculate the spike rate in Hz of self non-vocal behaviors cut into the
% same duartion chuncks as vocalizations (for fair comparisons 
SelfNVBehav_rate = cell(length(BehavType),1);
for bb=1:length(BehavType)
    IndNVBehavS = find(contains(Data.What, BehavType{bb}).*contains(Data.Who, 'self'));
    if ~isempty(IndNVBehavS)
        SelfNVBehav_rate{bb} = nan(length(IndVoc),1);
        [DurVoc,~] = sort(Data.Duration(IndVoc),'descend');
        DurNVB = Data.Duration(IndNVBehavS);
        UsedTime = zeros(length(IndNVBehavS),1);
        for oo=1:length(IndVoc)
            [GoodIdx,DurNVB,UsedTime] = recurDur(DurNVB, DurVoc(oo),UsedTime);
            if ~isempty(GoodIdx)
                SAT = Data.SpikesArrivalTimes_Behav{IndNVBehavS(GoodIdx)};
                NumSpikes = sum((SAT>=(UsedTime(GoodIdx)-DurVoc(oo))).*(SAT<(UsedTime(GoodIdx)+DurVoc(oo))));
                SelfNVBehav_rate{bb}(oo) = NumSpikes./DurVoc(oo).*10^3;
            else % Could not find an extract of non-vocal behavior long enough to extract data
                SelfNVBehav_rate{bb}(oo) = NaN;
            end
        end
        
    else
        SelfNVBehav_rate{bb} = NaN;
    end
end

% Calculate the spike rate in Hz of others non-vocal behaviors
OthersNVBehav_rate = cell(length(BehavType),1);
for bb=1:length(BehavType)
    IndNVBehavO = find(contains(Data.What, BehavType{bb}).*~contains(Data.Who, 'self'));
    if ~isempty(IndNVBehavO)
        OthersNVBehav_rate{bb} = nan(length(IndVoc),1);
        [DurVoc,~] = sort(Data.Duration(IndVoc),'descend');
        DurNVB = Data.Duration(IndNVBehavO);
        UsedTime = zeros(length(IndNVBehavO),1);
        for oo=1:length(IndVoc)
            [GoodIdx,DurNVB,UsedTime] = recurDur(DurNVB, DurVoc(oo),UsedTime);
            if ~isempty(GoodIdx)
                SAT = Data.SpikesArrivalTimes_Behav{IndNVBehavO(GoodIdx)};
                NumSpikes = sum((SAT>=(UsedTime(GoodIdx)-DurVoc(oo))).*(SAT<(UsedTime(GoodIdx)+DurVoc(oo))));
                OthersNVBehav_rate{bb}(oo) = NumSpikes./DurVoc(oo).*10^3;
            else % Could not find an extract of non-vocal behavior long enough to extract data
                OthersNVBehav_rate{bb}(oo) = NaN;
            end
        end
    else
       OthersNVBehav_rate{bb}= NaN;
    end
end

% save data
SpikeRate.SelfCall_rate = SelfCall_rate;
SpikeRate.SelfCall_exptype = SelfCall_exptype;
SpikeRate.OthersCall_rate = OthersCall_rate;
SpikeRate.OthersCall_exptype = OthersCall_exptype;
SpikeRate.SelfNVBehav_rate = SelfNVBehav_rate;
SpikeRate.OthersNVBehav_rate = OthersNVBehav_rate;
SpikeRate.BehavType = BehavType;
save(FullDataSetFile, 'SpikeRate', '-append')

%% INTERNAL FUNCTION
function [GoodIdx,DurNVB,UsedTime] = recurDur(DurNVB, DurVoci,UsedTime)
DiffDur = DurNVB-DurVoci;
DiffDurPos = find(DiffDur>=0);
if ~isempty(DiffDurPos)
    [~,Imin]=min(DiffDur(DiffDurPos));
    GoodIdx = DiffDurPos(Imin);
    UsedTime(GoodIdx) = UsedTime(GoodIdx)+DurVoci;
    DurNVB(DiffDurPos(Imin)) = DurNVB(DiffDurPos(Imin)) - DurVoci;
else
    GoodIdx = [];
end
end

end
