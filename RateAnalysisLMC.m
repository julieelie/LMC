addpath(genpath('/Users/elie/Documents/CODE/GitHub/SoundAnalysisBats'));
addpath(genpath('/Users/elie/Documents/CODE/GitHub/LoggerDataProcessing'));
addpath(genpath('/Users/elie/Documents/CODE/GitHub/LMC'));
addpath(genpath('/Users/elie/Documents/CODE/GitHub/GeneralCode'));
addpath(genpath('/Users/elie/Documents/CODE/GitHub/TimeVaryingAcSemSTRF'));
DatFig=0; %Set to 1 to see input data figures for each cell
OutFig = 1;%Set to 1 to see output data figures for each cell

%Filename = '59834_20190611_SSS_1-97.mat';
% Filename = '59834_20190610_SSS_1-130.mat';
%Path = '/Users/elie/Documents/LMCResults/';
HDPath = '/Volumes/JulieE8T/LMCResults/';
ServerPath = '/Volumes/server_home/users/JulieE/LMC/LMCResults';
% Path = '/Users/elie/Documents/ManipBats/LMC/ResultsFiles/';
% Path = '/Users/elie/Google Drive/BatmanData/';




%% List all cells
% AllFiles = dir(fullfile(HDPath,'*8*.mat')); % Cells from both Cooper and Hodor
% AllFiles = dir(fullfile(HDPath,'*7*.mat')); % Cells from Hank
AllFiles = [dir(fullfile(HDPath,'65701*.mat')); dir(fullfile(HDPath,'59834*.mat')); dir(fullfile(HDPath,'11689*.mat'))];

Files2run = zeros(length(AllFiles),1);
for ff=1:length(AllFiles)
    fprintf(1, 'Cell %d/%d\n', ff, length(AllFiles))
    Cell = load(fullfile(AllFiles(ff).folder, AllFiles(ff).name), 'What');
    if (length(strfind(AllFiles(ff).name, '_'))==3) && isfield(Cell, 'What') % Cells without the what field are files of another type or non-compiled cells probably because they are just noise
        Files2run(ff) = 1;
    end
    clear Cell
end
CellsPath = AllFiles(logical(Files2run));
NCells = length(CellsPath);

%% Check cell stability between 2 sessions
OpvsFr_Ks2stat = nan(NCells,1);
OpvsFr_Ks2statevenodd = nan(NCells,1);
OpvsFr_Ks2stat_boot = cell(NCells,1);
OpvsFr_Ks2pval_boot = nan(NCells,1);
NHO = 6;
NHF = 5;
NH = NHO + NHF; % presume total duration of experiments (6 hours operant max, 5 hours free session max)
OpvsFr_zs_KDE = nan(NCells,NH*60); % contains the per cell normalized KDE (/min)
OpvsFr_meandiff_zs_KDE = nan(NCells,1);% mean difference zscored KDE between Operant and Free for each cell
Nbootstrap = 10000;
for cc=1:NCells
    
    %% load data
    CellPath = fullfile(CellsPath(cc).folder,CellsPath(cc).name);
    fprintf(1, '*** Cell %s %d/%d ***\n', CellPath, cc, NCells)
    Cell = load(CellPath,'QualitySSU', 'FreeBehavSession', 'OperantSession');
    % decide if cell is stable enough to compare free and operant sessions
    LogicalOp=logical(((Cell.QualitySSU.TimePoints(2:end)-60/2)>0) .* ((Cell.QualitySSU.TimePoints(2:end)-60/2)<diff(Cell.OperantSession)));
    LogicalFr=logical(((Cell.QualitySSU.TimePoints(2:end)-60/2)>(Cell.FreeBehavSession(1)- Cell.OperantSession(1))) .* ((Cell.QualitySSU.TimePoints(2:end)-60/2)<(Cell.FreeBehavSession(2)- Cell.OperantSession(1))));
    IndicesOp = find(LogicalOp);
    IndicesFr = find(LogicalFr);
    if isnan(nanmean(Cell.QualitySSU.KDE))
        warning('Error with that cell, the time varying rate calculated by sanitary_check_perSSfile is Nan')
        continue
    end
    ZscoreKDE = (Cell.QualitySSU.KDE - nanmean(Cell.QualitySSU.KDE))/nanstd(Cell.QualitySSU.KDE);
    OpvsFr_zs_KDE(cc,(NHO*60-sum(LogicalOp)+1) : NHO*60) = ZscoreKDE(LogicalOp);
    OpvsFr_zs_KDE(cc,NHO*60+(1:sum(LogicalFr))) = ZscoreKDE(LogicalFr);
    OpvsFr_meandiff_zs_KDE(cc) = mean(ZscoreKDE(LogicalOp)) - mean(ZscoreKDE(LogicalFr));
    if isempty(IndicesFr)
        fprintf(1, 'No Free session for this cell, no stability to calculate!\n')
        continue
    end
    AllIndices = [IndicesOp IndicesFr];
    EvenIndOp = IndicesOp(rem(IndicesOp,2)==0);
    OddIndOp = IndicesOp(rem(IndicesOp,2)==1);
    EvenIndFr = IndicesFr(rem(IndicesFr,2)==0);
    OddIndFr = IndicesFr(rem(IndicesFr,2)==1);
    %     [H,KDE_OpvsFr_p(cc),CI,STATS] = ttest2(Cell.QualitySSU.KDE(logical(IndicesOp)), Cell.QualitySSU.KDE(logical(IndicesFr)));
    if sum(LogicalFr)&& sum(LogicalOp)
        [~,~,OpvsFr_Ks2stat(cc)] = kstest2(Cell.QualitySSU.KDE(LogicalOp), Cell.QualitySSU.KDE(LogicalFr));
        [~,~,OpvsFr_Ks2statevenodd(cc)] = kstest2(Cell.QualitySSU.KDE([EvenIndOp EvenIndFr]), Cell.QualitySSU.KDE([OddIndOp OddIndFr]));
        OpvsFr_Ks2stat_boot_local = cell(1,Nbootstrap);
        parfor bb=1:Nbootstrap
            Shuffled = AllIndices(randperm(length(AllIndices)));
            [~,~,OpvsFr_Ks2stat_boot_local{bb}] = kstest2(Cell.QualitySSU.KDE(Shuffled(1:length(IndicesOp))), Cell.QualitySSU.KDE(Shuffled(length(IndicesOp)+ (1:length(IndicesFr)))));
        end
        OpvsFr_Ks2stat_boot{cc} = [OpvsFr_Ks2stat_boot_local{:}];
    end
    OpvsFr_Ks2pval_boot(cc) = sum(OpvsFr_Ks2stat_boot{cc}>=OpvsFr_Ks2stat(cc))/Nbootstrap;
    
%     figure(15);clf;subplot(1,2,1); histogram(Cell.QualitySSU.KDE(LogicalOp), 'Normalization', 'probability'); hold on; histogram(Cell.QualitySSU.KDE(LogicalFr), 'Normalization', 'probability')
%     xlabel('Spike Rate in Hz (bin=1min)')
%     ylabel('probability')
%     legend({'operant' 'Free'})
%     title('Spike rate distribution')
%     subplot(1,2,2)
%     histogram(OpvsFr_Ks2stat_boot{cc}, 'FaceColor', 'k'); hold on; VObs = vline(OpvsFr_Ks2stat(cc), 'g-'); VObs.LineWidth=2; hold on; VEO = vline(OpvsFr_Ks2statevenodd(cc), 'r-'); VEO.LineWidth=2;
%     title(sprintf('KS difference distribution p=%.2f',OpvsFr_Ks2pval_boot(cc)))
%     pause(1)
end
save(fullfile(HDPath, 'KDE_SessionStability.mat'), 'CellsPath', 'OpvsFr_Ks2stat', 'OpvsFr_Ks2statevenodd', 'OpvsFr_Ks2stat_boot','OpvsFr_Ks2pval_boot', 'OpvsFr_meandiff_zs_KDE', 'OpvsFr_zs_KDE')
save(fullfile(ServerPath, 'KDE_SessionStability.mat'), 'CellsPath', 'OpvsFr_Ks2stat', 'OpvsFr_Ks2statevenodd', 'OpvsFr_Ks2stat_boot','OpvsFr_Ks2pval_boot', 'OpvsFr_meandiff_zs_KDE', 'OpvsFr_zs_KDE')
%%
NonNan = ~isnan(OpvsFr_Ks2pval_boot);
figure(13);
clf
subplot(1,3,1)
scatter(OpvsFr_Ks2stat(NonNan)-OpvsFr_Ks2statevenodd(NonNan), OpvsFr_Ks2stat(NonNan)-cellfun(@mean, OpvsFr_Ks2stat_boot(NonNan)),30,[OpvsFr_Ks2pval_boot(NonNan)>0.01 zeros(size(OpvsFr_Ks2pval_boot(NonNan))) zeros(size(OpvsFr_Ks2pval_boot(NonNan)))], 'filled'); 
xlabel('KS Obs - KS evenOdd'); ylabel('KS Obs - KS Bootstrap')
subplot(1,3,2)
LocalDat = OpvsFr_Ks2stat(NonNan)-cellfun(@mean, OpvsFr_Ks2stat_boot(NonNan));
histogram(LocalDat(OpvsFr_Ks2pval_boot(NonNan)<=0.01), 'BinWidth',0.05); hold on; histogram(LocalDat(OpvsFr_Ks2pval_boot(NonNan)>0.01), 'BinWidth',0.05);xlabel('KS stat Obs - Bootstraped'); ylabel('# cells'); legend('Affected','Unaffected')
title('Effect of session')
subplot(1,3,3)
LocalDat = OpvsFr_Ks2stat(NonNan)-OpvsFr_Ks2statevenodd(NonNan);
histogram(LocalDat(OpvsFr_Ks2pval_boot(NonNan)<=0.01), 'BinWidth',0.05); hold on; histogram(LocalDat(OpvsFr_Ks2pval_boot(NonNan)>0.01), 'BinWidth',0.05);xlabel('KS stat Obs - EvenOdd'); ylabel('# cells'); legend('Affected','Unaffected')
suplabel(sprintf('Cell change of rate with session (N=%d): %d unaffected cells', NCells, sum(OpvsFr_Ks2pval_boot(NonNan)>0.01)))

FigRate = figure(14);
clf
subplot(1,6,1:5)
[~,SortI] = sort(OpvsFr_meandiff_zs_KDE(NonNan)); %%% FIX THIS!!!!
OpvsFr_zs_KDE_NoNan = OpvsFr_zs_KDE(NonNan,:);
imagesc(OpvsFr_zs_KDE_NoNan(SortI,:))
colormap("parula")
colorbar
Clim = caxis;
AmpC = diff(Clim);
caxis(Clim(1) + [0 AmpC*0.75])
xlabel('Time (min)')
ylabel('Cell #')
title('Zscored-rate (/min)')
FigRate.CurrentAxes.XTick = [NHO*60/2 NHO*60+NHF*60/2];
FigRate.CurrentAxes.XTickLabel = {'Operant' 'Free'};
subplot(1,6,6)
OpvsFr_Ks2stat_NoNan = OpvsFr_Ks2stat(NonNan);
OpvsFr_Ks2stat_boot_NoNan = OpvsFr_Ks2stat_boot(NonNan);
OpvsFr_Ks2pval_boot_NoNan = OpvsFr_Ks2pval_boot(NonNan);
scatter(OpvsFr_Ks2stat_NoNan(SortI) - cellfun(@mean, OpvsFr_Ks2stat_boot_NoNan(SortI)), 1:length(OpvsFr_Ks2stat_NoNan), 20, [OpvsFr_Ks2pval_boot_NoNan(SortI)>0.01 zeros(size(OpvsFr_Ks2pval_boot_NoNan)) zeros(size(OpvsFr_Ks2pval_boot_NoNan))], 'filled')
set(gca, 'YDir', 'reverse')
xlabel('Cell Change of Rate Distribution (KS stat diff)')

%% Gather cell average rate around vocalizations
% Matrix of zscored average spike rate -/+5s vocalization onset Operant vocalizer,
% Operant Listner, Free session vocalizer, Free session Listner
Delay=4500;%(in ms)
ZSRate_Operant_Self = nan(NCells, 2*Delay+1);
ZSRate_Operant_Others = ZSRate_Operant_Self;
ZSRate_Free_Self = ZSRate_Operant_Self;
ZSRate_Free_Others = ZSRate_Operant_Self;
ZSRate_Operant_Self_Tr = ZSRate_Operant_Self;
ZSRate_Operant_Others_Tr = ZSRate_Operant_Self;
ZSRate_Free_Self_Tr = ZSRate_Operant_Self;
ZSRate_Free_Others_Tr = ZSRate_Operant_Self;
ZSRate_Operant_Self_NonTr = ZSRate_Operant_Self;
ZSRate_Operant_Others_NonTr = ZSRate_Operant_Self;
ZSRate_Free_Self_NonTr = ZSRate_Operant_Self;
ZSRate_Free_Others_NonTr = ZSRate_Operant_Self;
MaxPeakTime = nan(NCells,12);
MaxRate = nan(NCells,12);
CoherenceInfo = nan(NCells,12);
CoherenceInfoLow = nan(NCells,12);
CellsName = cell(NCells,1);
for cc=1:NCells
    % load data
    CellPath = fullfile(CellsPath(cc).folder,CellsPath(cc).name);
    fprintf(1, '*** Cell %s %d/%d ***\n', CellPath, cc, NCells)
    CellsName{cc} = CellPath;
    Cell = load(CellPath,'KDE_onset', 'MotorCoherenceOperant','MotorCoherenceFree', 'AuditoryCoherenceOperant', 'AuditoryCoherenceFree', 'MotorCoherenceOperantTrill', 'MotorCoherenceFreeTrill', 'AuditoryCoherenceOperantTrill', 'AuditoryCoherenceFreeTrill','MotorCoherenceOperantBark', 'MotorCoherenceFreeBark', 'AuditoryCoherenceOperantBark', 'AuditoryCoherenceFreeBark');
    if isfield(Cell.KDE_onset, 'SelfVocOp')
        Xt = selectXt(Cell.KDE_onset.SelfVocOp, Delay);
        ZSRate_Operant_Self(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.SelfVocOp(1,Xt));
        MaxRate(cc,1) = max(Cell.KDE_onset.SelfVocOp(1,Xt));
        [~ , MaxPeakTime(cc,1)] = max(ZSRate_Operant_Self(cc,:));
    end
    if isfield(Cell, 'MotorCoherenceOperant') && ~isfield(Cell.MotorCoherenceOperant, 'Error')
        CoherenceInfo(cc,1) = Cell.MotorCoherenceOperant.Info;
        CoherenceInfoLow(cc,1) = Cell.MotorCoherenceOperant.Info_low;
    end

    if isfield(Cell.KDE_onset, 'SelfVocFr')
        Xt = selectXt(Cell.KDE_onset.SelfVocFr, Delay);
        ZSRate_Free_Self(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.SelfVocFr(1,Xt));
        MaxRate(cc,2) = max(Cell.KDE_onset.SelfVocFr(1,Xt));
        [~ , MaxPeakTime(cc,2)] = max(ZSRate_Free_Self(cc,:));
    end
    if isfield(Cell, 'MotorCoherenceFree') && ~isfield(Cell.MotorCoherenceFree, 'Error')
        CoherenceInfo(cc,2) = Cell.MotorCoherenceFree.Info;
        CoherenceInfoLow(cc,2) = Cell.MotorCoherenceFree.Info_low;
    end

    if isfield(Cell.KDE_onset, 'OthersVocOp')
        Xt = selectXt(Cell.KDE_onset.OthersVocOp, Delay);
        ZSRate_Operant_Others(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.OthersVocOp(1,Xt));
        MaxRate(cc,3) = max(Cell.KDE_onset.OthersVocOp(1,Xt));
        [~ , MaxPeakTime(cc,3)] = max(ZSRate_Operant_Others(cc,:));
    end
    if isfield(Cell, 'AuditoryCoherenceOperant') && ~isfield(Cell.AuditoryCoherenceOperant, 'Error')
        CoherenceInfo(cc,3) = Cell.AuditoryCoherenceOperant.Info;
        CoherenceInfoLow(cc,3) = Cell.AuditoryCoherenceOperant.Info_low;
    end

    if isfield(Cell.KDE_onset, 'OthersVocFr')
        Xt = selectXt(Cell.KDE_onset.OthersVocFr, Delay);
        ZSRate_Free_Others(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.OthersVocFr(1,Xt));
        MaxRate(cc,4) = max(Cell.KDE_onset.OthersVocFr(1,Xt));
        [~ , MaxPeakTime(cc,4)] = max(ZSRate_Free_Others(cc,:));
    end
    if isfield(Cell, 'AuditoryCoherenceFree') && ~isfield(Cell.AuditoryCoherenceFree, 'Error')
        CoherenceInfo(cc,4) = Cell.AuditoryCoherenceFree.Info;
        CoherenceInfoLow(cc,4) = Cell.AuditoryCoherenceFree.Info_low;
    end



    if isfield(Cell.KDE_onset, 'SelfTrOp')
        Xt = selectXt(Cell.KDE_onset.SelfTrOp, Delay);
        ZSRate_Operant_Self_Tr(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.SelfTrOp(1,Xt));
        MaxRate(cc,5) = max(Cell.KDE_onset.SelfTrOp(1,Xt));
        [~ , MaxPeakTime(cc,5)] = max(ZSRate_Operant_Self_Tr(cc,:));
    end
    if isfield(Cell, 'MotorCoherenceOperantTrill') && ~isfield(Cell.MotorCoherenceOperantTrill, 'Error')
        CoherenceInfo(cc,5) = Cell.MotorCoherenceOperantTrill.Info;
        CoherenceInfoLow(cc,5) = Cell.MotorCoherenceOperantTrill.Info_low;
    end

    if isfield(Cell.KDE_onset, 'SelfTrFr')
        Xt = selectXt(Cell.KDE_onset.SelfTrFr, Delay);
        ZSRate_Free_Self_Tr(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.SelfTrFr(1,Xt));
        MaxRate(cc,6) = max(Cell.KDE_onset.SelfTrFr(1,Xt));
        [~ , MaxPeakTime(cc,6)] = max(ZSRate_Free_Self_Tr(cc,:));
    end
    if isfield(Cell, 'MotorCoherenceFreeTrill') && ~isfield(Cell.MotorCoherenceFreeTrill, 'Error')
        CoherenceInfo(cc,6) = Cell.MotorCoherenceFreeTrill.Info;
        CoherenceInfoLow(cc,6) = Cell.MotorCoherenceFreeTrill.Info_low;
    end

    if isfield(Cell.KDE_onset, 'OthersTrOp')
        Xt = selectXt(Cell.KDE_onset.OthersTrOp, Delay);
        ZSRate_Operant_Others_Tr(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.OthersTrOp(1,Xt));
        MaxRate(cc,7) = max(Cell.KDE_onset.OthersTrOp(1,Xt));
        [~ , MaxPeakTime(cc,7)] = max(ZSRate_Operant_Others_Tr(cc,:));
    end
    if isfield(Cell, 'AuditoryCoherenceOperantTrill') && ~isfield(Cell.AuditoryCoherenceOperantTrill, 'Error')
        CoherenceInfo(cc,7) = Cell.AuditoryCoherenceOperantTrill.Info;
        CoherenceInfoLow(cc,7) = Cell.AuditoryCoherenceOperantTrill.Info_low;
    end

    if isfield(Cell.KDE_onset, 'OthersTrFr')
        Xt = selectXt(Cell.KDE_onset.OthersTrFr, Delay);
        ZSRate_Free_Others_Tr(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.OthersTrFr(1,Xt));
        MaxRate(cc,8) = max(Cell.KDE_onset.OthersTrFr(1,Xt));
        [~ , MaxPeakTime(cc,8)] = max(ZSRate_Free_Others_Tr(cc,:));
    end
    if isfield(Cell, 'AuditoryCoherenceFreeTrill') && ~isfield(Cell.AuditoryCoherenceFreeTrill, 'Error')
        CoherenceInfo(cc,8) = Cell.AuditoryCoherenceFreeTrill.Info;
        CoherenceInfoLow(cc,8) = Cell.AuditoryCoherenceFreeTrill.Info_low;
    end


    if isfield(Cell.KDE_onset, 'SelfBaOp')
        Xt = selectXt(Cell.KDE_onset.SelfBaOp, Delay);
        ZSRate_Operant_Self_NonTr(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.SelfBaOp(1,Xt));
        MaxRate(cc,9) = max(Cell.KDE_onset.SelfBaOp(1,Xt));
        [~ , MaxPeakTime(cc,9)] = max(ZSRate_Operant_Self_NonTr(cc,:));
    end
    if isfield(Cell, 'MotorCoherenceOperantBark') && ~isfield(Cell.MotorCoherenceOperantBark, 'Error')
        CoherenceInfo(cc,9) = Cell.MotorCoherenceOperantBark.Info;
        CoherenceInfoLow(cc,9) = Cell.MotorCoherenceOperantBark.Info_low;
    end

    if isfield(Cell.KDE_onset, 'SelfBaFr')
        Xt = selectXt(Cell.KDE_onset.SelfBaFr, Delay);
        ZSRate_Free_Self_NonTr(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.SelfBaFr(1,Xt));
        MaxRate(cc,10) = max(Cell.KDE_onset.SelfBaFr(1,Xt));
        [~ , MaxPeakTime(cc,10)] = max(ZSRate_Free_Self_NonTr(cc,:));
    end
    if isfield(Cell, 'MotorCoherenceFreeBark') && ~isfield(Cell.MotorCoherenceFreeBark, 'Error')
        CoherenceInfo(cc,10) = Cell.MotorCoherenceFreeBark.Info;
        CoherenceInfoLow(cc,10) = Cell.MotorCoherenceFreeBark.Info_low;
    end

    if isfield(Cell.KDE_onset, 'OthersBaOp')
        Xt = selectXt(Cell.KDE_onset.OthersBaOp, Delay);
        ZSRate_Operant_Others_NonTr(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.OthersBaOp(1,Xt));
        MaxRate(cc,11) = max(Cell.KDE_onset.OthersBaOp(1,Xt));
        [~ , MaxPeakTime(cc,11)] = max(ZSRate_Operant_Others_NonTr(cc,:));
    end
    if isfield(Cell, 'AuditoryCoherenceOperantBark') && ~isfield(Cell.AuditoryCoherenceOperantBark, 'Error')
        CoherenceInfo(cc,11) = Cell.AuditoryCoherenceOperantBark.Info;
        CoherenceInfoLow(cc,11) = Cell.AuditoryCoherenceOperantBark.Info_low;
    end

    if isfield(Cell.KDE_onset, 'OthersBaFr')
        Xt = selectXt(Cell.KDE_onset.OthersBaFr, Delay);
        ZSRate_Free_Others_NonTr(cc,1:sum(Xt)) = myzscore(Cell.KDE_onset.OthersBaFr(1,Xt));
        MaxRate(cc,12) = max(Cell.KDE_onset.OthersBaFr(1,Xt));
        [~ , MaxPeakTime(cc,12)] = max(ZSRate_Free_Others_Tr(cc,:));
    end
    if isfield(Cell, 'AuditoryCoherenceFreeBark') && ~isfield(Cell.AuditoryCoherenceFreeBark, 'Error')
        CoherenceInfo(cc,12) = Cell.AuditoryCoherenceFreeBark.Info;
        CoherenceInfoLow(cc,12) = Cell.AuditoryCoherenceFreeBark.Info_low;
    end
end

 %% Plot the matrices and their averages
% Nnan = find(~isnan(MaxPeakTime(:,1)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,1)); % ordering the cells according to when they reach their max firing rate
% T=myplotmatrix(ZSRate_Operant_Self(Nnan(IndSort),:),Delay);
% title(T,'ZS Rate Self Operant')
% set(gca, 'CLim',[-5 9])

% Nnan = find(~isnan(MaxPeakTime(:,2)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,2)); % ordering the cells according to when they reach their max firing rate
% T=myplotmatrix(ZSRate_Free_Self(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Self Free')
% set(gca, 'CLim',[-5 9])
% 
% Nnan = find(~isnan(MaxPeakTime(:,3)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,3)); % ordering the cells according to when they reach their max firing rate
% T=myplotmatrix(ZSRate_Operant_Others(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Others Operant')
% set(gca, 'CLim',[-5 9])
% 
% Nnan = find(~isnan(MaxPeakTime(:,4)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,4)); % ordering the cells according to when they reach their max firing rate
% T=myplotmatrix(ZSRate_Free_Others(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Others Free')
% set(gca, 'CLim',[-5 9])
% 
% Nnan = find(~isnan(MaxPeakTime(:,5)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,5)); % ordering the cells according to when they reach their max firing rate
% [T, IM]=myplotmatrix(ZSRate_Operant_Self_Tr(Nnan(IndSort),:),Delay)
% title(T,'ZS Rate Self Operant Trills')
% set(IM, 'CLim',[-4 9])

Nnan = find(~isnan(MaxPeakTime(:,6)));
[~,IndSort] = sort(MaxPeakTime(Nnan,6)); % ordering the cells according to when they reach their max firing rate
[T,IM]=myplotmatrix(ZSRate_Free_Self_Tr(Nnan(IndSort),:), Delay)
title(T,'ZS Rate Self Free Trills')
set(IM, 'CLim',[-4 9])
% 
% Nnan = find(~isnan(MaxPeakTime(:,7)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,7)); % ordering the cells according to when they reach their max firing rate
% [T, IM]=myplotmatrix(ZSRate_Operant_Others_Tr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Others Operant Trills')
% set(IM, 'CLim',[-4 9])
% 
Nnan = find(~isnan(MaxPeakTime(:,8)));
[~,IndSort] = sort(MaxPeakTime(Nnan,8)); % ordering the cells according to when they reach their max firing rate
[T,IM]=myplotmatrix(ZSRate_Free_Others_Tr(Nnan(IndSort),:), Delay)
title(T,'ZS Rate Others Free Trills')
set(IM, 'CLim',[-4 9])

% Nnan = find(~isnan(MaxPeakTime(:,9)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,9)); % ordering the cells according to when they reach their max firing rate
% [T, IM]=myplotmatrix(ZSRate_Operant_Self_NonTr(Nnan(IndSort),:),Delay)
% title(T,'ZS Rate Self Operant Non-Trills')
% set(IM, 'CLim',[-4 9])

% Nnan = find(~isnan(MaxPeakTime(:,10)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,10)); % ordering the cells according to when they reach their max firing rate
% T=myplotmatrix(ZSRate_Free_Self_NonTr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Self Free Non-Trills')
% set(gca, 'CLim',[-5 9])
% % 
% Nnan = find(~isnan(MaxPeakTime(:,11)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,11)); % ordering the cells according to when they reach their max firing rate
% [T, IM]=myplotmatrix(ZSRate_Operant_Others_NonTr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Others Operant NonTrills')
% set(IM, 'CLim',[-4 9])
% 
% Nnan = find(~isnan(MaxPeakTime(:,12)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,12)); % ordering the cells according to when they reach their max firing rate
% T=myplotmatrix(ZSRate_Free_Others_NonTr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Others Free NonTrills')
% set(gca, 'CLim',[-5 9])

%% plot responses to trills from cells where the bats PRODUCE BOTH during
% Operant and Free session
Nnan = find((~isnan(MaxPeakTime(:,5))).*(~isnan(MaxPeakTime(:,6)))); % Cells for which we have both Trill production during Operant and Free sessions
[~,IndSort] = sort(MaxPeakTime(Nnan,5)); % ordering the cells according to when they reach their max firing rate in Operant Session
% T=myplotmatrix(ZSRate_Operant_Self_Tr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Self Operant Trills')
% set(gca, 'CLim',[-4 5.5])
% T=myplotmatrix(ZSRate_Free_Self_Tr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Self Free Trills')
% set(gca, 'CLim',[-4 5.5])
% [~,IndSort] = sort(MaxPeakTime(Nnan,6)); 
% T=myplotmatrix(ZSRate_Free_Self_Tr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Self Free Trills')
% set(gca, 'CLim',[-4 5.5])

% % Get the pearson correlation coefficient between session for the same
% % cells
% [Rho, PvalRho] = corr(ZSRate_Operant_Self_Tr(Nnan(IndSort),:)', ZSRate_Free_Self_Tr(Nnan(IndSort),:)');
% Rho = diag(Rho);
% PvalRho = diag(PvalRho);
% figure(18)
% clf
% [T, IM1, IM2, L]=myplotmatrix2(ZSRate_Operant_Self_Tr(Nnan(IndSort),:),Rho,Delay);
% % set(gca, 'XLim', [0 3])
% title(IM1,'ZS Rate')
% title(IM2, 'Rho')
% title(T, 'Operant Self Trills only')
% set(IM1, 'CLim',[-4 5.5])
% set(IM2, 'CLim',[0 1])
% set(L, 'YLim', [-0.5 0.5])
% 
% % histogram of correlation coefficients across the KDE -4.5sec ->4.5 sec
% figure(19)
% clf
% histogram(Rho,-1:0.05:1);
% hold on
% VL = vline(nanmean(Rho));
% VL.LineWidth=2;
% ylabel('# cells')
% xlabel('Rho')
% title('Correlation of time varying rates around Trill onset (4500ms) between Operant and Free')

% figures with correlation coefficients across the KDE -.5sec ->.5 sec
Delay2 = 500;
Timepoints = (Delay+(-Delay2:Delay2));
T=myplotmatrix(ZSRate_Operant_Self_Tr(Nnan(IndSort),Timepoints), Delay2)
title(T,'ZS Rate Self Operant Trills')
set(gca, 'CLim',[-4 5.5])
T=myplotmatrix(ZSRate_Free_Self_Tr(Nnan(IndSort),Timepoints), Delay2)
title(T,'ZS Rate Self Free Trills')
set(gca, 'CLim',[-4 5.5])

% [Rho, PvalRho] = corr(ZSRate_Operant_Self_Tr(Nnan(IndSort),Timepoints)', ZSRate_Free_Self_Tr(Nnan(IndSort),Timepoints)');
% Rho = diag(Rho);
% PvalRho = diag(PvalRho);
% figure(20)
% clf
% [T, IM1, IM2, L]=myplotmatrix2(ZSRate_Operant_Self_Tr(Nnan(IndSort),Timepoints),Rho,Delay2);
% % set(gca, 'XLim', [0 3])
% title(IM1,'ZS Rate')
% title(IM2, 'Rho')
% title(T, 'Operant Self Trills only')
% set(IM1, 'CLim',[-4 5.5])
% set(IM2, 'CLim',[0 1])
% set(L, 'YLim', [-0.5 0.5])
% 
% figure(21)
% clf
% histogram(Rho,-1:0.05:1);
% hold on
% VL = vline(nanmean(Rho));
% VL.LineWidth=2;
% ylabel('# cells')
% xlabel('Rho')
% title('Correlation of time varying rates around Trill onset (500ms) between Operant and Free')

% NS = histcounts(Rho(PvalRho>=0.01),-1:0.05:1)';
% Sig = histcounts(Rho(PvalRho<0.01),-1:0.05:1)';
% bar([NS Sig], 'stacked')

%% plot responses to trills from cells where the bats HEAR BOTH during
% Operant and Free session
% Nnan = find((~isnan(MaxPeakTime(:,7))).*(~isnan(MaxPeakTime(:,8)))); % Cells for which we have both Trill production during Operant and Free sessions
% [~,IndSort] = sort(MaxPeakTime(Nnan,7)); % ordering the cells according to when they reach their max firing rate in Operant Session
% T=myplotmatrix(ZSRate_Operant_Others_Tr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Others Operant Trills')
% set(gca, 'CLim',[-3 8])
% T=myplotmatrix(ZSRate_Free_Others_Tr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Others Free Trills')
% set(gca, 'CLim',[-3 8])
% [~,IndSort] = sort(MaxPeakTime(Nnan,8)); 
% T=myplotmatrix(ZSRate_Free_Others_Tr(Nnan(IndSort),:), Delay)
% title(T,'ZS Rate Others Free Trills')
% set(gca, 'CLim',[-3 8])

%% Matrices with barplot of Info on Coherence
Delay2 = 500;
Timepoints = (Delay+(-Delay2:Delay2));

% Nnan = find(~isnan(MaxPeakTime(:,1)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,1)); % ordering the cells according to when they reach their max firing rate
% figure(7)
% clf
% [T, IM1, IM2, L]=myplotmatrix2(ZSRate_Operant_Self(Nnan(IndSort),Timepoints),CoherenceInfo(Nnan(IndSort),1),Delay2);
% title(IM1,'ZS Rate')
% title(IM2,'Info (bit/s)')
% title(T, 'Self Operant')
% set(IM1, 'CLim',[-5 9])
% set(IM2, 'CLim',[0 3])
% set(L, 'YLim', [-0.6 0.8])
% 
% Nnan = find(~isnan(MaxPeakTime(:,2)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,2)); % ordering the cells according to when they reach their max firing rate
% figure(8)
% clf
% [T, IM1, IM2,L]=myplotmatrix2(ZSRate_Free_Self(Nnan(IndSort),Timepoints),CoherenceInfo(Nnan(IndSort),2),Delay2);
% % set(gca, 'XLim', [0 3])
% title(IM1,'ZS Rate')
% title(IM2, 'Info (bit/s)')
% title(T, 'Self Free')
% set(IM1, 'CLim',[-5 9])
% set(IM2, 'CLim',[0 3])
% set(L, 'YLim', [-0.05 0.6])
% 
% Nnan = find(~isnan(MaxPeakTime(:,3)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,3)); % ordering the cells according to when they reach their max firing rate
% figure(9)
% clf
% [T, IM1, IM2,L]=myplotmatrix2(ZSRate_Operant_Others(Nnan(IndSort),Timepoints),CoherenceInfo(Nnan(IndSort),3),Delay2);
% % set(gca, 'XLim', [0 3])
% title(IM1,'ZS Rate')
% title(IM2, 'Info (bit/s)')
% title(T, 'Operant Others')
% set(IM1, 'CLim',[-5 9])
% set(IM2, 'CLim',[0 3])
% set(L, 'YLim', [-0.6 0.8])
% 
% Nnan = find(~isnan(MaxPeakTime(:,4)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,4)); % ordering the cells according to when they reach their max firing rate
% figure(10)
% clf
% [T, IM1, IM2, L]=myplotmatrix2(ZSRate_Free_Others(Nnan(IndSort),Timepoints),CoherenceInfo(Nnan(IndSort),4),Delay2);
% % set(gca, 'XLim', [0 3])
% title(IM1,'ZS Rate')
% title(IM2, 'Info (bit/s)')
% title(T, 'Free Others')
% set(IM1, 'CLim',[-5 9])
% set(IM2, 'CLim',[0 3])
% set(L, 'YLim', [-0.05 0.6])

% Nnan = find(~isnan(MaxPeakTime(:,5)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,5)); % ordering the cells according to when they reach their max firing rate
% figure(11)
% clf
% [T, IM1, IM2, L]=myplotmatrix2(ZSRate_Operant_Self_Tr(Nnan(IndSort),Timepoints),CoherenceInfo(Nnan(IndSort),5),Delay2);
% % set(gca, 'XLim', [0 3])
% title(IM1,'ZS Rate')
% title(IM2, 'Info (bit/s)')
% title(T, 'Operant Self Trills only')
% set(IM1, 'CLim',[-5 9])
% set(IM2, 'CLim',[0 4])
% set(L, 'YLim', [-0.4 0.4])
% 
% Nnan = find(~isnan(MaxPeakTime(:,7)));
% [~,IndSort] = sort(MaxPeakTime(Nnan,7)); % ordering the cells according to when they reach their max firing rate
% figure(12)
% clf
% [T, IM1, IM2, L]=myplotmatrix2(ZSRate_Operant_Others_Tr(Nnan(IndSort),Timepoints),CoherenceInfo(Nnan(IndSort),7),Delay2);
% % set(gca, 'XLim', [0 3])
% title(IM1,'ZS Rate')
% title(IM2, 'Info (bit/s)')
% title(T, 'Operant Others Trills only')
% set(IM1, 'CLim',[-5 9])
% set(IM2, 'CLim',[0 4])
% set(L, 'YLim', [-0.4 0.4])

Nnan = find(~isnan(MaxPeakTime(:,9)));
[~,IndSort] = sort(MaxPeakTime(Nnan,9)); % ordering the cells according to when they reach their max firing rate
figure(13)
clf
[T, IM1, IM2, L]=myplotmatrix2(ZSRate_Operant_Self_NonTr(Nnan(IndSort),Timepoints),CoherenceInfo(Nnan(IndSort),9),Delay2);
% set(gca, 'XLim', [0 3])
title(IM1,'ZS Rate')
title(IM2, 'Info (bit/s)')
title(T, 'Operant Self non-Trills only')
set(IM1, 'CLim',[-3 8])
set(IM2, 'CLim',[0 6])
set(L, 'YLim', [-0.6 1.5])

Nnan = find(~isnan(MaxPeakTime(:,11)));
[~,IndSort] = sort(MaxPeakTime(Nnan,11)); % ordering the cells according to when they reach their max firing rate
figure(14)
clf
[T, IM1, IM2, L]=myplotmatrix2(ZSRate_Operant_Others_NonTr(Nnan(IndSort),Timepoints),CoherenceInfo(Nnan(IndSort),11),Delay2);
% set(gca, 'XLim', [0 3])
title(IM1,'ZS Rate')
title(IM2, 'Info (bit/s)')
title(T, 'Operant Others non-Trills only')
set(IM1, 'CLim',[-3 8])
set(IM2, 'CLim',[0 6])
set(L, 'YLim', [-0.6 1.5])


%% Distributions of MaxTime
figure()
histogram(MaxPeakTime(:,1))
hold on
histogram(MaxPeakTime(:,3))
legend('Operant Production', 'Operant Perception')

figure()
subplot(1,3,1:2)
    histogram(MaxPeakTime(:,1)-Delay,'BinWidth',10, 'FaceColor','r', 'normalization', 'cumcount')
    hold on
    histogram(MaxPeakTime(:,3)-Delay, 'BinWidth',10, 'FaceColor', 'b','normalization', 'cumcount')
    xlabel('Time of Peak spike rate (ms)')
    ylabel('Cumulative count of cells')
    hold on
    L = vline(0,'k-');
    L.LineWidth = 2;
    legend('Operant Production', 'Operant Perception')
    subplot(1,3,3)
    histogram(MaxPeakTime(:,1)-Delay,'BinWidth',10, 'FaceColor','r', 'normalization', 'cumcount')
    hold on
    histogram(MaxPeakTime(:,3)-Delay, 'BinWidth',10, 'FaceColor', 'b','normalization', 'cumcount')
    xlabel('Time of Peak spike rate (ms)')
    ylabel('Cumulative count of cells')
    hold on
    L = vline(0,'k-');
    L.LineWidth = 2;
    legend('Operant Production', 'Operant Perception')
    xlim([-500 500])
    ylim([0 350])

    

    Tbl =  table([Auditory.Info(logical(Row65701.*SigA)); Motor.Info(logical((Row59834 + Row11689).*SigM))], [repmat('A', sum(logical(Row65701.*SigA)),1);repmat('M', sum(logical((Row59834 + Row11689).*SigM)),1)], [repmat(65701,sum(logical(Row65701.*SigA)),1); Motor.BatID(logical((Row59834 + Row11689).*SigM))], 'VariableNames', {'Info', 'AorM','ID'});
    lme = fitlme(Tbl, 'Info~AorM+(1|ID)')
    lme2 = fitlme(Tbl, 'Info~1+(1|ID)')
    TestLME=compare(lme2, lme)
    if Trill==1
        title(sprintf('Vocalizing (n=%d in 2 bats) vs hearing (n= %d in 1 bat) Trills in %s (LME Info~AorM+(1|Subject) p=%.3f)', sum(logical((Row59834 + Row11689).*SigM)), sum(logical(Row65701.*SigA)), SessionChoice, TestLME.pValue(2)))
    elseif Trill==0
        title(sprintf('Vocalizing (n=%d in 2 bats) vs hearing (n= %d in 1 bat) any call in %s (LME Info~AorM+(1|Subject) p=%.3f)', sum(logical((Row59834 + Row11689).*SigM)), sum(logical(Row65701.*SigA)), SessionChoice, TestLME.pValue(2)))
    elseif Trill==2
        title(sprintf('Vocalizing (n=%d in 2 bats) vs hearing (n= %d in 1 bat) non-Trill in %s (LME Info~AorM+(1|Subject) p=%.3f)', sum(logical((Row59834 + Row11689).*SigM)), sum(logical(Row65701.*SigA)), SessionChoice, TestLME.pValue(2)))
    end
    legend({'Perception'; 'Production'})




function [T, IM1]=myplotmatrix(DAT1, Delay)
    figure()
    T=tiledlayout(5,1);
    nexttile([4,1])
    imagesc(DAT1)
    IM1 = gca;
    colorbar('northoutside')
    colormap('cool')
    hold on
    line([Delay+1 Delay+1], [0 size(DAT1,1)],'Color','y','LineStyle','--', 'LineWidth',2)
    hold on
    ylabel('Units')
    xlabel('Time (s)')
    xticks([1 Delay+1 2*Delay+1])
    xticklabels([-Delay 0 Delay])
    hold off
    nexttile
    shadedErrorBar(-Delay:Delay,mean(DAT1),std(DAT1)./(size(DAT1,1))^0.5,{'k-','LineWidth',2});
    ylabel('ZS spike rate')
    xlabel('Time (s)')
    xlim([-Delay Delay])
end

function [T, IM1, IM2, L]=myplotmatrix2(DAT1,DAT2, Delay)
    T=tiledlayout(5,5);
    nexttile([4,4])
    imagesc(DAT1);
    IM1 = gca;
    colorbar(IM1,'northoutside')
    colormap(IM1,'cool')
    hold on
    line([Delay+1 Delay+1], [0 size(DAT1,1)],'Color','y','LineStyle','--', 'LineWidth',2)
    hold on
    ylabel('Units')
    xlabel('Time (s)')
    xticks([1 Delay+1 2*Delay+1])
    xticklabels([-Delay 0 Delay])
    hold off
    nexttile(21,[1 4])
    shadedErrorBar(-Delay:Delay,mean(DAT1),std(DAT1)./(size(DAT1,1))^0.5,{'k-','LineWidth',2});
    ylabel('ZS spike rate')
    xlabel('Time (s)')
    xlim([-Delay Delay])
    L = gca;
    nexttile(5, [4 1])
    imagesc(DAT2)
    IM2 = gca;
    colormap(IM2,"bone")
    colorbar(IM2,'Eastoutside')
%     barh(DAT2, 0.7, 'k');
%     set(gca, 'YDir', 'reverse')
end


function DATOUT = myzscore(DAT)
    DATOUT = round(DAT-mean(DAT),5)./std(DAT);
end

function Xt = selectXt(Dat, Delay)
    Xt = logical((Dat(2,:)>=-Delay) .* (Dat(2,:)<=Delay));
end