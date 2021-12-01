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
AllFiles = dir(fullfile(HDPath,'*7*.mat')); % Cells from both Hank
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
OpvsFr_zs_KDE = nan(NCells,600); % contains the per cell normalized KDE (/min)
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
    ZscoreKDE = (Cell.QualitySSU.KDE - nanmean(Cell.QualitySSU.KDE))/nanstd(Cell.QualitySSU.KDE);
    OpvsFr_zs_KDE(cc,(300-sum(LogicalOp)+1) : 300) = ZscoreKDE(LogicalOp);
    OpvsFr_zs_KDE(cc,300+(1:sum(LogicalFr))) = ZscoreKDE(LogicalFr);
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
figure();
subplot(1,3,1)
scatter(OpvsFr_Ks2stat-OpvsFr_Ks2statevenodd, OpvsFr_Ks2stat-cellfun(@mean, OpvsFr_Ks2stat_boot),30,[OpvsFr_Ks2pval_boot>0.01 zeros(size(OpvsFr_Ks2pval_boot)) zeros(size(OpvsFr_Ks2pval_boot))], 'filled'); 
xlabel('KS Obs - KS evenOdd'); ylabel('KS Obs - KS Bootstrap')
subplot(1,3,2)
LocalDat = OpvsFr_Ks2stat-cellfun(@mean, OpvsFr_Ks2stat_boot);
histogram(LocalDat(OpvsFr_Ks2pval_boot<=0.01), 'BinWidth',0.05); hold on; histogram(LocalDat(OpvsFr_Ks2pval_boot>0.01), 'BinWidth',0.05);xlabel('KS stat Obs - Bootstraped'); ylabel('# cells'); legend('Unstable','Stable')
subplot(1,3,3)
LocalDat = OpvsFr_Ks2stat-OpvsFr_Ks2statevenodd;
histogram(LocalDat(OpvsFr_Ks2pval_boot<=0.01), 'BinWidth',0.05); hold on; histogram(LocalDat(OpvsFr_Ks2pval_boot>0.01), 'BinWidth',0.05);xlabel('KS stat Obs - EvenOdd'); ylabel('# cells'); legend('Unstable','Stable')
suplabel(sprintf('Cell stability (N=%d): %d stable cells', NCells, sum(OpvsFr_Ks2pval_boot>0.01)))

FigRate = figure();
subplot(1,6,1:5)
[~,SortI] = sort(OpvsFr_meandiff_zs_KDE);
imagesc(OpvsFr_zs_KDE(SortI,1:500))
colormap("parula")
colorbar
Clim = caxis;
AmpC = diff(Clim);
caxis(Clim(1) + [0 AmpC*0.75])
xlabel('Time (min)')
ylabel('Cell #')
title('zscored-rate')
FigRate.CurrentAxes.XTick = [150 400];
FigRate.CurrentAxes.XTickLabel = {'Operant' 'Free'};
subplot(1,6,6)
scatter(OpvsFr_Ks2stat(SortI) - cellfun(@mean, OpvsFr_Ks2stat_boot(SortI)), 1:length(OpvsFr_Ks2stat), 20, [OpvsFr_Ks2pval_boot(SortI)>0.01 zeros(size(OpvsFr_Ks2pval_boot)) zeros(size(OpvsFr_Ks2pval_boot))], 'filled')
set(gca, 'YDir', 'reverse')
xlabel('Cell Instability (KS stat diff)')

%% Plot cell average rate around vocalizations
% Matrix of zscored average spike rate -/+5s vocalization onset Operant vocalizer,
% Operant Listner, Free session vocalizer, Free session Listner
Delay=4500;%(in ms)
ZSRate_Operant_Self = nan(NCells, 2*Delay+1);
ZSRate_Operant_Others = ZSRate_Operant_Self;
ZSRate_Free_Self = ZSRate_Operant_Self;
ZSRate_Free_Others = ZSRate_Operant_Self;
MaxPeakTime = nan(NCells,4);
for cc=1:NCells
    %% load data
    CellPath = fullfile(CellsPath(cc).folder,CellsPath(cc).name);
    fprintf(1, '*** Cell %s %d/%d ***\n', CellPath, cc, NCells)
    Cell = load(CellPath,'KDE_onset');
    if isfield(KDE_onset, 'SelfVocOp')
        Xt = selectXt(Cell.KDE_onset.SelfVocOp, Delay);
        ZSRate_Operant_Self(cc,1:length(Xt)) = zscore(Cell.KDE_onset.SelfVocOp(1,Xt));
        [~ , MaxPeakTime(cc,1)] = max(ZSRate_Operant_Self(cc,:));
    end
    if isfield(KDE_onset, 'SelfVocFr')
        Xt = selectXt(Cell.KDE_onset.SelfVocFr, Delay);
        ZSRate_Free_Self(cc,1:length(Xt)) = zscore(Cell.KDE_onset.SelfVocFr(1,Xt));
        [~ , MaxPeakTime(cc,2)] = max(ZSRate_Free_Self(cc,:));
    end
    if isfield(KDE_onset, 'OthersVocOp')
        Xt = selectXt(Cell.KDE_onset.OthersVocOp, Delay);
        ZSRate_Operant_Others(cc,1:length(Xt)) = zscore(Cell.KDE_onset.OthersVocOp(1,Xt));
        [~ , MaxPeakTime(cc,3)] = max(ZSRate_Operant_Others(cc,:));
    end
    if isfield(KDE_onset, 'OthersVocFr')
        Xt = selectXt(Cell.KDE_onset.OthersVocFr, Delay);
        ZSRate_Free_Others(cc,1:length(Xt)) = zscore(Cell.KDE_onset.OthersVocFr(1,Xt));
        [~ , MaxPeakTime(cc,4)] = max(ZSRate_Free_Others(cc,:));
    end
end

figure()
subplot(5,1,1:4)
Nnan = find(~isnan(MaxPeakTime(:,1)));
[~,IndSort] = sort(MaxPeakTime(Nnan,1)); % ordering the cells according to when they reach their max firing rate
imagesc(ZSRate_Operant_Self(Nnan(IndSort),:))
colorbar()
colormap('cool')
hold on
line([Delay+1 Delay+1], [0 length(Nnan)],'Color','y','LineStyle','--', 'LineWidth',2)
hold on
ylabel('Units')
xlabel('Time (s)')
xticks([1 Delay+1 2*Delay+1])
xticklabels([-Delay 0 Delay])
title('ZS Rate Self Operant')
hold off
subplot(5,1,5)
shadedErrorBar([],ZSRate_Operant_Self(Nnan(IndSort),:),{@mean,@std},{'k-','LineWidth',2});
ylabel('ZS spike rate')
xlabel('Time (s)')
xticks([1 Delay+1 2*Delay+1])
xticklabels([-Delay 0 Delay])


figure()
subplot(5,1,1:4)
Nnan = find(~isnan(MaxPeakTime(:,2)));
[~,IndSort] = sort(MaxPeakTime(Nnan,2)); % ordering the cells according to when they reach their max firing rate
imagesc(ZSRate_Free_Self(Nnan(IndSort),:))
colorbar()
colormap('cool')
hold on
line([Delay+1 Delay+1], [0 length(Nnan)],'Color','y','LineStyle','--', 'LineWidth',2)
hold on
ylabel('Units')
xlabel('Time (s)')
xticks([1 Delay+1 2*Delay+1])
xticklabels([-Delay 0 Delay])
title('ZS Rate Self Free')
hold off
subplot(5,1,5)
shadedErrorBar([],ZSRate_Free_Self(Nnan(IndSort),:),{@mean,@std},{'k-','LineWidth',2});
ylabel('ZS spike rate')
xlabel('Time (s)')
xticks([1 Delay+1 2*Delay+1])
xticklabels([-Delay 0 Delay])

figure()
subplot(5,1,1:4)
Nnan = find(~isnan(MaxPeakTime(:,3)));
[~,IndSort] = sort(MaxPeakTime(Nnan,3)); % ordering the cells according to when they reach their max firing rate
imagesc(ZSRate_Operant_Others(Nnan(IndSort),:))
colorbar()
colormap('cool')
hold on
line([Delay+1 Delay+1], [0 length(Nnan)],'Color','y','LineStyle','--', 'LineWidth',2)
hold on
ylabel('Units')
xlabel('Time (s)')
xticks([1 Delay+1 2*Delay+1])
xticklabels([-Delay 0 Delay])
title('ZS Rate Others Operant')
hold off
subplot(5,1,5)
shadedErrorBar([],ZSRate_Operant_Others(Nnan(IndSort),:),{@mean,@std},{'k-','LineWidth',2});
ylabel('ZS spike rate')
xlabel('Time (s)')
xticks([1 Delay+1 2*Delay+1])
xticklabels([-Delay 0 Delay])


figure()
subplot(5,1,1:4)
Nnan = find(~isnan(MaxPeakTime(:,4)));
[~,IndSort] = sort(MaxPeakTime(Nnan,4)); % ordering the cells according to when they reach their max firing rate
imagesc(ZSRate_Free_Others(Nnan(IndSort),:))
colorbar()
colormap('cool')
hold on
line([Delay+1 Delay+1], [0 length(Nnan)],'Color','y','LineStyle','--', 'LineWidth',2)
hold on
ylabel('Units')
xlabel('Time (s)')
xticks([1 Delay+1 2*Delay+1])
xticklabels([-Delay 0 Delay])
title('ZS Rate Others Free')
hold off
subplot(5,1,5)
shadedErrorBar([],ZSRate_Free_Others(Nnan(IndSort),:),{@mean,@std},{'k-','LineWidth',2});
ylabel('ZS spike rate')
xlabel('Time (s)')
xticks([1 Delay+1 2*Delay+1])
xticklabels([-Delay 0 Delay])

function Xt = selectXt(Dat, Delay)
    Xt = logical((Dat(2,:)>=-Delay) .* (Dat(2,:)<=Delay));
end