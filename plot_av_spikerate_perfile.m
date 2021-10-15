function plot_av_spikerate_perfile(InputDataFile, OutputPath, Fun)
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


%% Plot the figure
ScatterMarkerSz = 30;
MeanMarkerSize = 14;
ColorBR = [0 0.4470 0.7410];
ColorSelf = 'k';
ColorOthers = [0.9290 0.6940 0.1250];
Fig=figure(); 
Nevents = zeros(8,1);
% Plot the spike rate in Hz of self calls Operant
Ind = contains(SpikeRate.SelfCall_exptype, 'O');
Nevents(1:2) = sum(Ind).*ones(2,1);
swarmchart(ones(sum(Ind),1),Fun(SpikeRate.SelfCall_rate(Ind,1)),ScatterMarkerSz,ColorSelf,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
hold on
errorbar(1,mean(Fun(SpikeRate.SelfCall_rate(Ind,1))),std(Fun(SpikeRate.SelfCall_rate(Ind,1)))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
hold on
swarmchart(2.*ones(sum(Ind),1),Fun(SpikeRate.SelfCall_rate(Ind,2)),ScatterMarkerSz,ColorBR,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
hold on
errorbar(2,mean(Fun(SpikeRate.SelfCall_rate(Ind,2))),std(Fun(SpikeRate.SelfCall_rate(Ind,2)))/sum(Ind)^0.5, 'dc','MarkerSize',MeanMarkerSize,'MarkerFaceColor','c')
hold on


% Plot the spike rate in Hz of self calls Free session
Ind = contains(SpikeRate.SelfCall_exptype, 'F');
Nevents(3:4) = sum(Ind).*ones(2,1);
swarmchart(3.*ones(sum(Ind),1),Fun(SpikeRate.SelfCall_rate(Ind,1)),ScatterMarkerSz,ColorSelf,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
hold on
errorbar(3,mean(Fun(SpikeRate.SelfCall_rate(Ind,1))),std(Fun(SpikeRate.SelfCall_rate(Ind,1)))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
hold on
swarmchart(4.*ones(sum(Ind),1),Fun(SpikeRate.SelfCall_rate(Ind,2)),ScatterMarkerSz,ColorBR,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
hold on
errorbar(4,nanmean(Fun(SpikeRate.SelfCall_rate(Ind,2))),nanstd(Fun(SpikeRate.SelfCall_rate(Ind,2)))/sum(Ind)^0.5, 'dc','MarkerSize',MeanMarkerSize,'MarkerFaceColor','c')
hold on


% Plot the spike rate in Hz of others calls Operant
if ~isempty(SpikeRate.OthersCall_exptype)
    % This is useless as the companion at rarely vocalize
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
    
    
    
    % Plot the spike rate in Hz of others calls Free session
    Ind = contains(SpikeRate.OthersCall_exptype, 'F');
    Nevents(5:6) = sum(Ind).*ones(2,1);
    swarmchart(5.*ones(sum(Ind),1),Fun(SpikeRate.OthersCall_rate(Ind,1)),ScatterMarkerSz,ColorOthers,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    hold on
    errorbar(5,mean(Fun(SpikeRate.OthersCall_rate(Ind,1))),std(Fun(SpikeRate.OthersCall_rate(Ind,1)))/sum(Ind)^0.5, 'sr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
    hold on
    swarmchart(6.*ones(sum(Ind),1),Fun(SpikeRate.OthersCall_rate(Ind,2)),ScatterMarkerSz,ColorBR,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    hold on
    errorbar(6,nanmean(Fun(SpikeRate.OthersCall_rate(Ind,2))),nanstd(Fun(SpikeRate.OthersCall_rate(Ind,2)))/sum(Ind)^0.5, 'sc','MarkerSize',MeanMarkerSize,'MarkerFaceColor','c')
    hold on
    
end
LegendVoc = {'Self-Voc-Operant' 'bSelf-Voc-Operant' 'Self-Voc-Free' 'bSelf-Voc-Free' 'Others-Voc-Free' 'bOthers-Voc-Free'};
% Indicate the number of events
text(1,-0.5,sprintf('%d',Nevents(1)))
text(2,-0.5,sprintf('%d',Nevents(2)))
text(3,-0.5,sprintf('%d',Nevents(3)))
text(4,-0.5,sprintf('%d',Nevents(4)))
text(5,-0.5,sprintf('%d',Nevents(5)))
text(6,-0.5,sprintf('%d',Nevents(6)))
% text(7,-0.5,sprintf('%d',Nevents(7)))
% text(8,-0.5,sprintf('%d',Nevents(8)))

if length(SpikeRate.BehavType)==0
    keyboard
end

% Plot Self non vocal behavior (Free session)
LegendSNVB = cell(1,length(SpikeRate.BehavType));
for bb=1:length(SpikeRate.BehavType)
    if ~isnan(SpikeRate.SelfNVBehav_rate{bb})
%         Jitter_local = rand([size(SpikeRate.SelfNVBehav_rate{bb},1),1]).*0.4-0.2;
        swarmchart((bb+6).*ones(size(SpikeRate.SelfNVBehav_rate{bb},1),1),Fun(SpikeRate.SelfNVBehav_rate{bb}),ScatterMarkerSz,ColorSelf,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
        hold on
        errorbar(bb+6,mean(Fun(SpikeRate.SelfNVBehav_rate{bb})),std(Fun(SpikeRate.SelfNVBehav_rate{bb}))/size(SpikeRate.SelfNVBehav_rate{bb},1)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
        hold on
        
        if contains(SpikeRate.BehavType{bb}, 'teeth')
            LegendSNVB{bb} = sprintf('Self-%s','nailbiting');
        else
            LegendSNVB{bb} = sprintf('Self-%s',SpikeRate.BehavType{bb});
        end
        text(bb+6,-0.5,sprintf('%d',size(SpikeRate.SelfNVBehav_rate{bb},1)))
    else
        if contains(SpikeRate.BehavType{bb}, 'teeth')
            LegendSNVB{bb} = sprintf('Self-%s','nailbiting');
        else
            LegendSNVB{bb} = sprintf('Self-%s',SpikeRate.BehavType{bb});
        end
        text(bb+6,-0.5,'0')
    end
end

% Plot Others non vocal behavior (Free session)
LegendONVB = cell(1,length(SpikeRate.BehavType));
for bb=1:length(SpikeRate.BehavType)
    if ~isnan(SpikeRate.OthersNVBehav_rate{bb})
%         Jitter_local = rand([size(SpikeRate.OthersNVBehav_rate{bb},1),1]).*0.4-0.2;
        swarmchart((bb+6+length(SpikeRate.BehavType)).*ones(size(SpikeRate.OthersNVBehav_rate{bb},1),1),Fun(SpikeRate.OthersNVBehav_rate{bb}),ScatterMarkerSz,ColorOthers,'o','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
        hold on
        errorbar(bb+6+length(SpikeRate.BehavType),mean(Fun(SpikeRate.OthersNVBehav_rate{bb})),std(Fun(SpikeRate.SelfNVBehav_rate{bb}))/size(SpikeRate.SelfNVBehav_rate{bb},1)^0.5, 'sr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
        hold on
        
        if contains(SpikeRate.BehavType{bb}, 'teeth')
            LegendONVB{bb} = sprintf('Others-%s','nailbiting');
        else
            LegendONVB{bb} = sprintf('Others-%s',SpikeRate.BehavType{bb});
        end
        text(bb+6+length(SpikeRate.BehavType),-0.5,sprintf('%d',size(SpikeRate.OthersNVBehav_rate{bb},1)))
    else
        if contains(SpikeRate.BehavType{bb}, 'teeth')
            LegendONVB{bb} = sprintf('Others-%s','nailbiting');
        else
            LegendONVB{bb} = sprintf('Others-%s',SpikeRate.BehavType{bb});
        end
        text(bb+6+length(SpikeRate.BehavType),-0.5,'0')
    end
        
end

ylabel('Spike rate (Hz)')
FullLegend = [LegendVoc LegendSNVB LegendONVB];
set(gca,'XTick', 1:length(FullLegend),'XTickLabel',FullLegend)
title(sprintf('Mean rate\n%s on %s Tetrode %s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}))
xlim(gca,[0 6+2*length(SpikeRate.BehavType)+1])

orient(Fig,'landscape')
Fig.PaperPositionMode = 'auto';
set(Fig,'PaperOrientation','landscape');
Fig.Units = 'inches';
Fig.Children.YLim(1)=-1;
Fig.Children.YLim(2) = min(100, max(1,Fig.Children.YLim(2)));
FigPosition = Fig.Position;
Fig.Position = [FigPosition(1:2) 10 4];
xtickangle(20)
print(Fig,fullfile(OutputPath,sprintf('%s_%s_SS%s_%s-%s_MeanRateScatter.pdf', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2})),'-dpdf','-fillpage')
close all
end
