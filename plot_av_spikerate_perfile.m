function [Stats]=plot_av_spikerate_perfile(InputDataFile, OutputPath, Fun)
if nargin<3
    Fun = @(x)x;
end
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
% Get the SS Quality
NeuralInputID{3} = DataFile(strfind(DataFile, '_SS')+3);

% Get the subject ID
SubjectID = DataFile(1:5);

% Input
FullDataSetFile = fullfile(OutputPath, sprintf('%s_%s_SS%s_%s-%s.mat', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2}));
load(FullDataSetFile, 'SpikeRate');

% Output
Stats = table('Size',[10 6],'VariableTypes',["string", 'double', 'double','double','int16', 'int16'],'VariableNames',{'Test', 'p-value','t-stat','DF','n1','n2'});
TestCount = 0;
%% Plot the figure
ScatterMarkerSz = 30;
MeanMarkerSize = 14;
ColorBR = [0 0.4470 0.7410];
ColorSelf = 'k';
ColorOthers = [0.9290 0.6940 0.1250];
Fig=figure();
Nevents = zeros(8,1);
LegendVoc = {};
NaxisInd = -1;
% Plot the spike rate in Hz of self calls Operant
Ind = contains(SpikeRate.SelfCall_exptype, 'O');
if sum(Ind)>5
    NaxisInd = NaxisInd+2;
    Nevents(NaxisInd:NaxisInd+1) = sum(Ind).*ones(2,1);
    swarmchart(NaxisInd.*ones(sum(Ind),1),Fun(SpikeRate.SelfCall_rate(Ind,1)),ScatterMarkerSz,ColorSelf,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    hold on
    errorbar(NaxisInd,mean(Fun(SpikeRate.SelfCall_rate(Ind,1))),std(Fun(SpikeRate.SelfCall_rate(Ind,1)))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
    hold on
    swarmchart((NaxisInd+1).*ones(sum(Ind),1),Fun(SpikeRate.SelfCall_rate(Ind,2)),ScatterMarkerSz,ColorBR,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    hold on
    errorbar((NaxisInd+1),nanmean(Fun(SpikeRate.SelfCall_rate(Ind,2))),nanstd(Fun(SpikeRate.SelfCall_rate(Ind,2)))/sum(Ind)^0.5, 'dc','MarkerSize',MeanMarkerSize,'MarkerFaceColor','c')
    hold on
    [h,pSO,ci,stats] = ttest(SpikeRate.SelfCall_rate(Ind,1),SpikeRate.SelfCall_rate(Ind,2));
    fprintf(1,'Operant: Vocalization production vs background rate (1s 5s before call onset): p = %.4f   t=%.2f   DF = %.1f\n', pSO, stats.tstat, stats.df)
    LegendVoc = [LegendVoc{:} {'Self-Voc-Operant' 'bSelf-Voc-Operant'}];
    TestCount = TestCount+1;
    Stats(TestCount,:) = {'Voc-Operant-SelfVsBgd', pSO,stats.tstat, stats.df,sum(Ind),sum(Ind)};
end


% Plot the spike rate in Hz of self calls Free session
Ind = contains(SpikeRate.SelfCall_exptype, 'F');
if sum(Ind)>5
    NaxisInd = NaxisInd+2;
    Nevents(NaxisInd:(NaxisInd+1)) = sum(Ind).*ones(2,1);
    swarmchart(NaxisInd.*ones(sum(Ind),1),Fun(SpikeRate.SelfCall_rate(Ind,1)),ScatterMarkerSz,ColorSelf,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    hold on
    errorbar(NaxisInd,mean(Fun(SpikeRate.SelfCall_rate(Ind,1))),std(Fun(SpikeRate.SelfCall_rate(Ind,1)))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
    hold on
    swarmchart((NaxisInd+1).*ones(sum(Ind),1),Fun(SpikeRate.SelfCall_rate(Ind,2)),ScatterMarkerSz,ColorBR,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    hold on
    errorbar((NaxisInd+1),nanmean(Fun(SpikeRate.SelfCall_rate(Ind,2))),nanstd(Fun(SpikeRate.SelfCall_rate(Ind,2)))/sum(Ind)^0.5, 'dc','MarkerSize',MeanMarkerSize,'MarkerFaceColor','c')
    hold on
    [h,pSF,ci,stats] = ttest(SpikeRate.SelfCall_rate(Ind,1),SpikeRate.SelfCall_rate(Ind,2));
    fprintf(1,'Free: Vocalization production vs background rate (1s 5s before call onset): p = %.4f   t=%.2f   DF = %.1f\n', pSF, stats.tstat, stats.df)
    IndSelf = Ind;
    LegendVoc = [LegendVoc{:} {'Self-Voc-Free' 'bSelf-Voc-Free'}];
    TestCount = TestCount+1;
    Stats(TestCount,:) = {'Voc-Free-SelfVsBgd', pSF,stats.tstat, stats.df,sum(Ind),sum(Ind)};
end
% Plot the spike rate in Hz of others calls Operant
if ~isempty(SpikeRate.OthersCall_exptype)
    % This is a companion bat that hears but do not vocalize in Operant
    % conditioning
    %     Ind = contains(SpikeRate.OthersCall_exptype, 'O');
    %     Nevents(5:6) = sum(Ind).*ones(2,1);
    %     scatter(5.*ones(sum(Ind),1)+rand([sum(Ind),1]).*0.4-0.2,SpikeRate.OthersCall_rate(Ind,1),ScatterMarkerSz,'k','o','filled')
    %     hold on
    %     errorbar(5,mean(SpikeRate.OthersCall_rate(Ind,1)),std(SpikeRate.OthersCall_rate(Ind,1))/sum(Ind)^0.5, 'sr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
    %     hold on
    %     scatter(6.*ones(sum(Ind),1)+rand([sum(Ind),1]).*0.4-0.2,SpikeRate.OthersCall_rate(Ind,2),ScatterMarkerSz,'k','o','filled')
    %     hold on
    %     errorbar(6,mean(SpikeRate.OthersCall_rate(Ind,2)),std(SpikeRate.OthersCall_rate(Ind,2))/sum(Ind)^0.5, 'sc','MarkerSize',MeanMarkerSize,'MarkerFaceColor','c')
    %     hold on

    Ind = contains(SpikeRate.OthersCall_exptype, 'O');
    if sum(Ind)>5
        NaxisInd = NaxisInd+2;
        Nevents(NaxisInd:(NaxisInd+1)) = sum(Ind).*ones(2,1);
        swarmchart(NaxisInd.*ones(sum(Ind),1),Fun(SpikeRate.OthersCall_rate(Ind,1)),ScatterMarkerSz,ColorOthers,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
        hold on
        errorbar(NaxisInd,mean(Fun(SpikeRate.OthersCall_rate(Ind,1))),std(Fun(SpikeRate.OthersCall_rate(Ind,1)))/sum(Ind)^0.5, 'sr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
        hold on
        swarmchart((NaxisInd+1).*ones(sum(Ind),1),Fun(SpikeRate.OthersCall_rate(Ind,2)),ScatterMarkerSz,ColorBR,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
        hold on
        errorbar((NaxisInd+1),nanmean(Fun(SpikeRate.OthersCall_rate(Ind,2))),nanstd(Fun(SpikeRate.OthersCall_rate(Ind,2)))/sum(Ind)^0.5, 'sc','MarkerSize',MeanMarkerSize,'MarkerFaceColor','c')
        hold on
        [h,pOO,ci,stats] = ttest(SpikeRate.OthersCall_rate(Ind,1),SpikeRate.OthersCall_rate(Ind,2));
        fprintf(1,'Operant: Vocalization Others vs background rate (1s 5s before call onset): p = %.4f   t=%.2f   DF = %.1f\n', pOO, stats.tstat, stats.df)
        LegendVoc = [LegendVoc{:} {'Others-Voc-Operant' 'bOthers-Voc-Operant'}];
        TestCount = TestCount+1;
        Stats(TestCount,:) = {'Voc-Operant-OthersVsBgd', pOO,stats.tstat, stats.df,sum(Ind),sum(Ind)};
    end



    % Plot the spike rate in Hz of others calls Free session
    Ind = contains(SpikeRate.OthersCall_exptype, 'F');
    if sum(Ind)>5
        NaxisInd = NaxisInd +2;
        Nevents(NaxisInd:(NaxisInd+1)) = sum(Ind).*ones(2,1);
        swarmchart(NaxisInd.*ones(sum(Ind),1),Fun(SpikeRate.OthersCall_rate(Ind,1)),ScatterMarkerSz,ColorOthers,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
        hold on
        errorbar(NaxisInd,mean(Fun(SpikeRate.OthersCall_rate(Ind,1))),std(Fun(SpikeRate.OthersCall_rate(Ind,1)))/sum(Ind)^0.5, 'sr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
        hold on
        swarmchart((NaxisInd+1).*ones(sum(Ind),1),Fun(SpikeRate.OthersCall_rate(Ind,2)),ScatterMarkerSz,ColorBR,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
        hold on
        errorbar((NaxisInd+1),nanmean(Fun(SpikeRate.OthersCall_rate(Ind,2))),nanstd(Fun(SpikeRate.OthersCall_rate(Ind,2)))/sum(Ind)^0.5, 'sc','MarkerSize',MeanMarkerSize,'MarkerFaceColor','c')
        hold on
        [h,pOF,ci,stats] = ttest(SpikeRate.OthersCall_rate(Ind,1),SpikeRate.OthersCall_rate(Ind,2));
        fprintf(1,'Free: Vocalization Others vs background rate (1s 5s before call onset): p = %.4f   t=%.2f   DF = %.1f\n', pOF, stats.tstat, stats.df)
        TestCount = TestCount+1;
        Stats(TestCount,:) = {'Voc-Free-OthersVsBgd', pOF,stats.tstat, stats.df,sum(Ind),sum(Ind)};
        
        if exist('IndSelf', 'var') % minimum of 5 calls to get into the plot and t-test
            [h,p3,ci,stats] = ttest2(SpikeRate.SelfCall_rate(IndSelf,1), SpikeRate.OthersCall_rate(Ind,1));
            fprintf(1,'Free: Vocalization Self vs Others: p = %.4f   t=%.2f   DF = %.1f\n', p3, stats.tstat, stats.df)
            TestCount = TestCount+1;
            Stats(TestCount,:) = {'Voc-Free-SelfVsOthers', p3,stats.tstat, stats.df,sum(IndSelf),sum(Ind)};
        end
        LegendVoc = [LegendVoc{:} {'Others-Voc-Free' 'bOthers-Voc-Free'}];
    end


end

% Indicate the number of events
NaxisCount = NaxisInd+1;
for ee=1:(NaxisCount)
    text(ee,-0.5,sprintf('%d',Nevents(ee)))
end
%     text(2,-0.5,sprintf('%d',Nevents(2)))
%     text(3,-0.5,sprintf('%d',Nevents(3)))
%     text(4,-0.5,sprintf('%d',Nevents(4)))
% text(5,-0.5,sprintf('%d',Nevents(5)))
% text(6,-0.5,sprintf('%d',Nevents(6)))
% text(7,-0.5,sprintf('%d',Nevents(7)))
% text(8,-0.5,sprintf('%d',Nevents(8)))

if length(SpikeRate.BehavType)==0
    ylabel('Spike rate (Hz)')
    set(gca,'XTick', 1:length(LegendVoc),'XTickLabel',LegendVoc)
    title(sprintf('Mean rate\n%s on %s Tetrode %s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}))
    xlim(gca,[0 length(LegendVoc)+1])

    orient(Fig,'landscape')
    Fig.PaperPositionMode = 'auto';
    set(Fig,'PaperOrientation','landscape');
    Fig.Units = 'inches';
    Fig.Children.YLim(1)=-1;
    Fig.Children.YLim(2) = min(100, max(1,Fig.Children.YLim(2)));
    FigPosition = Fig.Position;
    Fig.Position = [FigPosition(1:2) 10 4];
    xtickangle(20)

    Ylim = get(gca, 'YLim');
    XInd = find(strcmp(LegendVoc, 'Self-Voc-Operant'));
    if ~isempty(XInd)
        line([XInd XInd+1], Ylim(2).*ones(2,1).*0.9, 'Color','k', 'LineWidth',2)
        if pSO<0.001
            text(XInd+0.2, Ylim(2).*0.9, '***', 'FontSize',25)
        elseif pSO<0.01
            text(XInd+0.3, Ylim(2).*0.9, '**', 'FontSize',25)
        elseif pSO<0.05
            text(XInd+0.4, Ylim(2).*0.9, '*', 'FontSize',25)
        else
            text(XInd+0.3, Ylim(2)*0.925, 'NS', 'FontSize',12)
        end
    end

    XInd = find(strcmp(LegendVoc, 'Self-Voc-Free'));
    if ~isempty(XInd)
        line([XInd XInd+1], Ylim(2).*ones(2,1).*0.9, 'Color','k', 'LineWidth',2)
        if pSF<0.001
            text(XInd+0.2, Ylim(2).*0.9, '***', 'FontSize',25)
        elseif pSF<0.01
            text(XInd+0.3, Ylim(2).*0.9, '**', 'FontSize',25)
        elseif pSF<0.05
            text(XInd+0.4, Ylim(2).*0.9, '*', 'FontSize',25)
        else
            text(XInd+0.3, Ylim(2)*0.925, 'NS', 'FontSize',12)
        end
    end

    XInd = find(strcmp(LegendVoc, 'Others-Voc-Operant'));
    if ~isempty(XInd)
        line([XInd XInd+1], Ylim(2).*ones(2,1).*0.9, 'Color','k', 'LineWidth',2)
        if pOO<0.001
            text(XInd+0.2, Ylim(2).*0.9, '***', 'FontSize',25)
        elseif pOO<0.01
            text(XInd+0.3, Ylim(2).*0.9, '**', 'FontSize',25)
        elseif pOO<0.05
            text(XInd+0.4, Ylim(2).*0.9, '*', 'FontSize',25)
        else
            text(XInd+0.3, Ylim(2)*0.925, 'NS', 'FontSize',12)
        end
    end

    XInd = find(strcmp(LegendVoc, 'Others-Voc-Free'));
    if ~isempty(XInd)
        line([XInd XInd+1], Ylim(2).*ones(2,1).*0.9, 'Color','k', 'LineWidth',2)
        if pOF<0.001
            text(XInd+0.2, Ylim(2).*0.9, '***', 'FontSize',25)
        elseif pOF<0.01
            text(XInd+0.3, Ylim(2).*0.9, '**', 'FontSize',25)
        elseif pOF<0.05
            text(XInd+0.4, Ylim(2).*0.9, '*', 'FontSize',25)
        else
            text(XInd+0.3, Ylim(2)*0.925, 'NS', 'FontSize',12)
        end

        XInd2 = find(strcmp(LegendVoc, 'Self-Voc-Free'));
        if ~isempty(XInd2)
            line([XInd2 XInd], Ylim(2).*ones(2,1).*0.95,'Color', 'k', 'LineWidth',2)
            if p3<0.001
                text(XInd2+0.7, Ylim(2).*0.95, '***', 'FontSize',25)
            elseif p3<0.01
                text(XInd2+0.8, Ylim(2).*0.95, '**', 'FontSize',25)
            elseif p3<0.05
                text(XInd2+0.9, Ylim(2).*0.95, '*', 'FontSize',25)
            else
                text(XInd2+0.8, Ylim(2)*0.975, 'NS', 'FontSize',12)
            end
        end
    end

%     keyboard
    print(Fig,fullfile(OutputPath,sprintf('%s_%s_SS%s_%s-%s_MeanRateScatter.pdf', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2})),'-dpdf','-fillpage')
    close all

else

    % Plot Self non vocal behavior (Free session)
    LegendSNVB = cell(1,length(SpikeRate.BehavType));
    for bb=1:length(SpikeRate.BehavType)
        if ~isnan(SpikeRate.SelfNVBehav_rate{bb})
            %         Jitter_local = rand([size(SpikeRate.SelfNVBehav_rate{bb},1),1]).*0.4-0.2;
            swarmchart((bb+NaxisCount).*ones(size(SpikeRate.SelfNVBehav_rate{bb},1),1),Fun(SpikeRate.SelfNVBehav_rate{bb}),ScatterMarkerSz,ColorSelf,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
            hold on
            errorbar(bb+NaxisCount,mean(Fun(SpikeRate.SelfNVBehav_rate{bb})),std(Fun(SpikeRate.SelfNVBehav_rate{bb}))/size(SpikeRate.SelfNVBehav_rate{bb},1)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
            hold on

            if contains(SpikeRate.BehavType{bb}, 'teeth')
                LegendSNVB{bb} = sprintf('Self-%s','nailbiting');
            else
                LegendSNVB{bb} = sprintf('Self-%s',SpikeRate.BehavType{bb});
            end
            text(bb+NaxisCount,-0.5,sprintf('%d',size(SpikeRate.SelfNVBehav_rate{bb},1)))
        else
            if contains(SpikeRate.BehavType{bb}, 'teeth')
                LegendSNVB{bb} = sprintf('Self-%s','nailbiting');
            else
                LegendSNVB{bb} = sprintf('Self-%s',SpikeRate.BehavType{bb});
            end
            text(bb+NaxisCount,-0.5,'0')
        end
    end

    % Plot Others non vocal behavior (Free session)
    LegendONVB = cell(1,length(SpikeRate.BehavType));
    for bb=1:length(SpikeRate.BehavType)
        if ~isnan(SpikeRate.OthersNVBehav_rate{bb})
            %         Jitter_local = rand([size(SpikeRate.OthersNVBehav_rate{bb},1),1]).*0.4-0.2;
            swarmchart((bb+NaxisCount+length(SpikeRate.BehavType)).*ones(size(SpikeRate.OthersNVBehav_rate{bb},1),1),Fun(SpikeRate.OthersNVBehav_rate{bb}),ScatterMarkerSz,ColorOthers,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
            hold on
            errorbar(bb+NaxisCount+length(SpikeRate.BehavType),mean(Fun(SpikeRate.OthersNVBehav_rate{bb})),std(Fun(SpikeRate.SelfNVBehav_rate{bb}))/size(SpikeRate.SelfNVBehav_rate{bb},1)^0.5, 'sr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
            hold on

            if contains(SpikeRate.BehavType{bb}, 'teeth')
                LegendONVB{bb} = sprintf('Others-%s','nailbiting');
            else
                LegendONVB{bb} = sprintf('Others-%s',SpikeRate.BehavType{bb});
            end
            text(bb+NaxisCount+length(SpikeRate.BehavType),-0.5,sprintf('%d',size(SpikeRate.OthersNVBehav_rate{bb},1)))
        else
            if contains(SpikeRate.BehavType{bb}, 'teeth')
                LegendONVB{bb} = sprintf('Others-%s','nailbiting');
            else
                LegendONVB{bb} = sprintf('Others-%s',SpikeRate.BehavType{bb});
            end
            text(bb+NaxisCount+length(SpikeRate.BehavType),-0.5,'0')
        end

    end

    ylabel('Spike rate (Hz)')
    FullLegend = [LegendVoc LegendSNVB LegendONVB];
    set(gca,'XTick', 1:length(FullLegend),'XTickLabel',FullLegend)
    title(sprintf('Mean rate\n%s on %s Tetrode %s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}))
    xlim(gca,[0 length(LegendVoc)+2*length(SpikeRate.BehavType)+1])

    orient(Fig,'landscape')
    Fig.PaperPositionMode = 'auto';
    set(Fig,'PaperOrientation','landscape');
    Fig.Units = 'inches';
    Fig.Children.YLim(1)=-1;
    Fig.Children.YLim(2) = min(100, max(1,Fig.Children.YLim(2)));
    FigPosition = Fig.Position;
    Fig.Position = [FigPosition(1:2) 10 4];
    xtickangle(20)


    Ylim = get(gca, 'YLim');
    XInd = find(strcmp(LegendVoc, 'Self-Voc-Operant'));
    if ~isempty(XInd)
        line([XInd XInd+1], Ylim(2).*ones(2,1).*0.9, 'Color','k', 'LineWidth',2)
        if pSO<0.001
            text(XInd+0.2, Ylim(2).*0.9, '***', 'FontSize',25)
        elseif pSO<0.01
            text(XInd+0.3, Ylim(2).*0.9, '**', 'FontSize',25)
        elseif pSO<0.05
            text(XInd+0.4, Ylim(2).*0.9, '*', 'FontSize',25)
        else
            text(XInd+0.3, Ylim(2)*0.925, 'NS', 'FontSize',12)
        end
    end

    
    XInd = find(strcmp(LegendVoc, 'Self-Voc-Free'));
    if ~isempty(XInd)
        line([XInd XInd+1], Ylim(2).*ones(2,1).*0.9, 'Color','k', 'LineWidth',2)
        if pSF<0.001
            text(XInd+0.2, Ylim(2).*0.9, '***', 'FontSize',25)
        elseif pSF<0.01
            text(XInd+0.3, Ylim(2).*0.9, '**', 'FontSize',25)
        elseif pSF<0.05
            text(XInd+0.4, Ylim(2).*0.9, '*', 'FontSize',25)
        else
            text(XInd+0.3, Ylim(2)*0.925, 'NS', 'FontSize',12)
        end
    end

    XInd = find(strcmp(LegendVoc, 'Others-Voc-Operant'));
    if ~isempty(XInd)
        line([XInd XInd+1], Ylim(2).*ones(2,1).*0.9, 'Color','k', 'LineWidth',2)
        if pOO<0.001
            text(XInd+0.2, Ylim(2).*0.9, '***', 'FontSize',25)
        elseif pOO<0.01
            text(XInd+0.3, Ylim(2).*0.9, '**', 'FontSize',25)
        elseif pOO<0.05
            text(XInd+0.4, Ylim(2).*0.9, '*', 'FontSize',25)
        else
            text(XInd+0.3, Ylim(2)*0.925, 'NS', 'FontSize',12)
        end
    end

    XInd = find(strcmp(LegendVoc, 'Others-Voc-Free'));
    if ~isempty(XInd)
        line([XInd XInd+1], Ylim(2).*ones(2,1).*0.9, 'Color','k', 'LineWidth',2)
        if pOF<0.001
            text(XInd+0.2, Ylim(2).*0.9, '***', 'FontSize',25)
        elseif pOF<0.01
            text(XInd+0.3, Ylim(2).*0.9, '**', 'FontSize',25)
        elseif pOF<0.05
            text(XInd+0.4, Ylim(2).*0.9, '*', 'FontSize',25)
        else
            text(XInd+0.3, Ylim(2)*0.925, 'NS', 'FontSize',12)
        end


        XInd2 = find(strcmp(LegendVoc, 'Self-Voc-Free'));
        if ~isempty(XInd2)
            line([XInd2 XInd], Ylim(2).*ones(2,1).*0.95,'Color', 'k', 'LineWidth',2)
            if p3<0.001
                text(XInd2+0.7, Ylim(2).*0.95, '***', 'FontSize',25)
            elseif p3<0.01
                text(XInd2+0.8, Ylim(2).*0.95, '**', 'FontSize',25)
            elseif p3<0.05
                text(XInd2+0.9, Ylim(2).*0.95, '*', 'FontSize',25)
            else
                text(XInd2+0.8, Ylim(2)*0.975, 'NS', 'FontSize',12)
            end
        end
    end

    % t-test Self chewing and self licking vs self quiet
    bbQuiet = find(contains(SpikeRate.BehavType, 'quiet'));
    bbChewing = find(contains(SpikeRate.BehavType, 'chewing'));
    bbLicking = find(contains(SpikeRate.BehavType, 'licking'));
    if ~isempty(bbQuiet) && ~isempty(bbChewing)
        [~, p4, ~,stats] = ttest2(SpikeRate.SelfNVBehav_rate{bbQuiet}, SpikeRate.SelfNVBehav_rate{bbChewing});
        fprintf(1,'Free: Chewing vs Quiet: p = %.4f   t=%.2f   DF = %.1f\n', p4, stats.tstat, stats.df)
        TestCount = TestCount+1;
        Stats(TestCount,:) = {'Self-Free-ChewingVsQuiet', p4,stats.tstat, stats.df,length(SpikeRate.SelfNVBehav_rate{bbChewing}),length(SpikeRate.SelfNVBehav_rate{bbQuiet})};
    end
    if ~isempty(bbQuiet) && ~isempty(bbLicking)
        [~, p5, ~,stats] = ttest2(SpikeRate.SelfNVBehav_rate{bbQuiet}, SpikeRate.SelfNVBehav_rate{bbLicking});
        fprintf(1,'Free: Licking vs Quiet: p = %.4f   t=%.2f   DF = %.1f\n', p5, stats.tstat, stats.df)
        TestCount = TestCount+1;
        Stats(TestCount,:) = {'Self-Free-LickingVsQuiet', p5,stats.tstat, stats.df,length(SpikeRate.SelfNVBehav_rate{bbLicking}),length(SpikeRate.SelfNVBehav_rate{bbQuiet})};
    end
    if ~isempty(bbQuiet) && exist('IndSelf', 'var')
        [~, p6, ~,stats] = ttest2(SpikeRate.SelfNVBehav_rate{bbQuiet}, SpikeRate.SelfCall_rate(IndSelf,1));
        fprintf(1,'Free: Vocalizing vs Quiet: p = %.4f   t=%.2f   DF = %.1f\n', p6, stats.tstat, stats.df)
        TestCount = TestCount+1;
        Stats(TestCount,:) = {'Self-Free-VocalizingVsQuiet', p6,stats.tstat, stats.df,sum(IndSelf),length(SpikeRate.SelfNVBehav_rate{bbQuiet})};
    end
    if ~isempty(bbQuiet) && ~isempty(bbChewing)
        line(length(LegendVoc)+[bbQuiet bbChewing], Ylim(2).*ones(2,1).*0.9, 'Color','k', 'LineWidth',2)
        TextX = length(LegendVoc)+min(bbQuiet,bbChewing)+abs(bbQuiet-bbChewing)/2;
        if p4<0.001
            text(TextX-0.3, Ylim(2).*0.9, '***', 'FontSize',25)
        elseif p4<0.01
            text(TextX-0.2, Ylim(2).*0.9, '**', 'FontSize',25)
        elseif p4<0.05
            text(TextX-0.1, Ylim(2).*0.9, '*', 'FontSize',25)
        else
            text(TextX-0.2, Ylim(2)*0.925, 'NS', 'FontSize',12)
        end
    end
    if ~isempty(bbQuiet) && ~isempty(bbLicking)
        line(length(LegendVoc)+[bbQuiet bbLicking], Ylim(2).*ones(2,1).*0.95, 'Color','k', 'LineWidth',2)
        TextX = length(LegendVoc)+min(bbQuiet, bbLicking) + abs(bbQuiet-bbLicking)/2;
        if p5<0.001
            text(TextX-0.3, Ylim(2).*0.95, '***', 'FontSize',25)
        elseif p5<0.01
            text(TextX-0.2, Ylim(2).*0.95, '**', 'FontSize',25)
        elseif p5<0.05
            text(TextX-0.1, Ylim(2).*0.95, '*', 'FontSize',25)
        else
            text(TextX-0.2, Ylim(2)*0.975, 'NS', 'FontSize',12)
        end
    end
    if ~isempty(bbQuiet) && exist('IndSelf', 'var')
        line([find(strcmp(LegendVoc, 'Self-Voc-Free')) length(LegendVoc)+bbQuiet], Ylim(2).*ones(2,1).*1,'Color', 'k', 'LineWidth',2)
        TextX = find(strcmp(LegendVoc, 'Self-Voc-Free')) + (3 + bbQuiet)/2;
        if p6<0.001
            text(TextX-0.3, Ylim(2).*1, '***', 'FontSize',25)
        elseif p6<0.01
            text(TextX-0.2, Ylim(2).*1, '**', 'FontSize',25)
        elseif p6<0.05
            text(TextX-0.1, Ylim(2).*1, '*', 'FontSize',25)
        else
            text(TextX-0.2, Ylim(2)*1.025, 'NS', 'FontSize',12)
        end
        Fig.Children.YLim(2) = Ylim(2)*1.05;
    end
    print(Fig,fullfile(OutputPath,sprintf('%s_%s_SS%s_%s-%s_MeanRateScatter.pdf', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2})),'-dpdf','-fillpage')
    close all
end
Stats = Stats(1:TestCount,:);
SpikeRate.Stats = Stats;
save(FullDataSetFile, 'SpikeRate', '-append');
end
