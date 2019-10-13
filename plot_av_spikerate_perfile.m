function plot_av_spikerate_perfile(InputDataFile, OutputPath)
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

% Input
FullDataSetFile = fullfile(OutputPath, sprintf('%s_%s_%s_SSU%s-%s.mat', SubjectID, Date,NeuralInputID{1},NeuralInputID{2}));
load(FullDataSetFile, 'SpikeRate');


%% Plot the figure
ScatterMarkerSz = 30;
MeanMarkerSize = 14;
Fig=figure();
Nevents = zeros(8,1);
% Plot the spike rate in Hz of self calls Operant
Ind = contains(SpikeRate.SelfCall_exptype, 'O');
Nevents(1:2) = sum(Ind).*ones(2,1);
errorbar(1,mean(SpikeRate.SelfCall_rate(Ind,1)),std(SpikeRate.SelfCall_rate(Ind,1))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
hold on
scatter(ones(sum(Ind),1),SpikeRate.SelfCall_rate(Ind,1),ScatterMarkerSz,'k','o','filled')
hold on
errorbar(2,mean(SpikeRate.SelfCall_rate(Ind,2)),std(SpikeRate.SelfCall_rate(Ind,2))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
hold on
scatter(2.*ones(sum(Ind),1),SpikeRate.SelfCall_rate(Ind,2),ScatterMarkerSz,'k','o','filled')
hold on

% Plot the spike rate in Hz of self calls Free session
Ind = contains(SpikeRate.SelfCall_exptype, 'F');
Nevents(3:4) = sum(Ind).*ones(2,1);
errorbar(3,mean(SpikeRate.SelfCall_rate(Ind,1)),std(SpikeRate.SelfCall_rate(Ind,1))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
hold on
scatter(3.*ones(sum(Ind),1),SpikeRate.SelfCall_rate(Ind,1),ScatterMarkerSz,'k','o','filled')
hold on
errorbar(4,mean(SpikeRate.SelfCall_rate(Ind,2)),std(SpikeRate.SelfCall_rate(Ind,2))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
hold on
scatter(4.*ones(sum(Ind),1),SpikeRate.SelfCall_rate(Ind,2),ScatterMarkerSz,'k','o','filled')
hold on

% Plot the spike rate in Hz of others calls Operant
if ~isempty(SpikeRate.OthersCall_exptype)
    Ind = contains(SpikeRate.OthersCall_exptype, 'O');
    Nevents(5:6) = sum(Ind).*ones(2,1);
    errorbar(5,mean(SpikeRate.OthersCall_rate(Ind,1)),std(SpikeRate.OthersCall_rate(Ind,1))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
    hold on
    scatter(5.*ones(sum(Ind),1),SpikeRate.OthersCall_rate(Ind,1),ScatterMarkerSz,'k','o','filled')
    hold on
    errorbar(6,mean(SpikeRate.OthersCall_rate(Ind,2)),std(SpikeRate.OthersCall_rate(Ind,2))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
    hold on
    scatter(6.*ones(sum(Ind),1),SpikeRate.OthersCall_rate(Ind,2),ScatterMarkerSz,'k','o','filled')
    hold on
    
    
    % Plot the spike rate in Hz of others calls Free session
    Ind = contains(SpikeRate.OthersCall_exptype, 'F');
    Nevents(7:8) = sum(Ind).*ones(2,1);
    errorbar(7,mean(SpikeRate.OthersCall_rate(Ind,1)),std(SpikeRate.OthersCall_rate(Ind,1))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
    hold on
    scatter(7.*ones(sum(Ind),1),SpikeRate.OthersCall_rate(Ind,1),ScatterMarkerSz,'k','o','filled')
    hold on
    errorbar(8,mean(SpikeRate.OthersCall_rate(Ind,2)),std(SpikeRate.OthersCall_rate(Ind,2))/sum(Ind)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
    hold on
    scatter(8.*ones(sum(Ind),1),SpikeRate.OthersCall_rate(Ind,2),ScatterMarkerSz,'k','o','filled')
    hold on
end
LegendVoc = {sprintf('S-Voc-O(%d)',Nevents(1)) sprintf('BS-Voc-O(%d)',Nevents(2)) sprintf('S-Voc-F(%d)',Nevents(3)) sprintf('BS-Voc-F(%d)',Nevents(4)) sprintf('O-Voc-O(%d)',Nevents(5)) sprintf('BO-Voc-O(%d)',Nevents(6)) sprintf('O-Voc-F(%d)',Nevents(7)) sprintf('BO-Voc-F(%d)',Nevents(8))};

LegendSNVB = cell(1,length(SpikeRate.BehavType));
for bb=1:length(SpikeRate.BehavType)
    if ~isnan(SpikeRate.SelfNVBehav_rate{bb})
        errorbar(bb+8,mean(SpikeRate.SelfNVBehav_rate{bb}),std(SpikeRate.SelfNVBehav_rate{bb})/size(SpikeRate.SelfNVBehav_rate{bb},1)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
        hold on
        scatter((bb+8).*ones(size(SpikeRate.SelfNVBehav_rate{bb},1),1),SpikeRate.SelfNVBehav_rate{bb},ScatterMarkerSz,'k','o','filled')
        hold on
        if contains(SpikeRate.BehavType{bb}, 'teeth')
            LegendSNVB{bb} = sprintf('S-%s(%d)','nailbiting',size(SpikeRate.SelfNVBehav_rate{bb},1));
        else
            LegendSNVB{bb} = sprintf('S-%s(%d)',SpikeRate.BehavType{bb},size(SpikeRate.SelfNVBehav_rate{bb},1));
        end
    else
        if contains(SpikeRate.BehavType{bb}, 'teeth')
            LegendSNVB{bb} = sprintf('S-%s(0)','nailbiting');
        else
            LegendSNVB{bb} = sprintf('S-%s(0)');
        end
    end
end

LegendONVB = cell(1,length(SpikeRate.BehavType));
for bb=1:length(SpikeRate.BehavType)
    if ~isnan(SpikeRate.OthersNVBehav_rate{bb})
        errorbar(bb+8+length(SpikeRate.BehavType),mean(SpikeRate.OthersNVBehav_rate{bb}),std(SpikeRate.SelfNVBehav_rate{bb})/size(SpikeRate.SelfNVBehav_rate{bb},1)^0.5, 'dr','MarkerSize',MeanMarkerSize,'MarkerFaceColor','r')
        hold on
        scatter((bb+8+length(SpikeRate.BehavType)).*ones(size(SpikeRate.OthersNVBehav_rate{bb},1),1),SpikeRate.OthersNVBehav_rate{bb},ScatterMarkerSz,'k','o','filled')
        hold on
        if contains(SpikeRate.BehavType{bb}, 'teeth')
            LegendONVB{bb} = sprintf('O-%s(%d)','nailbiting',size(SpikeRate.SelfNVBehav_rate{bb},1));
        else
            LegendONVB{bb} = sprintf('O-%s(%d)',SpikeRate.BehavType{bb},size(SpikeRate.SelfNVBehav_rate{bb},1));
        end
    else
        if contains(SpikeRate.BehavType{bb}, 'teeth')
            LegendONVB{bb} = sprintf('O-%s(0)','nailbiting');
        else
            LegendONVB{bb} = sprintf('O-%s(0)',SpikeRate.BehavType{bb});
        end
    end
        
end

ylabel('Spike rate (Hz)')
FullLegend = [LegendVoc LegendSNVB LegendONVB];
set(gca,'XTick', 1:length(FullLegend),'XTickLabel',FullLegend)
title(sprintf('Mean rate\n%s on %s Tetrode %s SU %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{2}))
xlim(gca,[0 8+2*length(SpikeRate.BehavType)+1])

orient(Fig,'landscape')
Fig.PaperPositionMode = 'auto';
set(Fig,'PaperOrientation','landscape');
print(Fig,fullfile(OutputPath,sprintf('%s_%s_%s_SSU%s-%s_MeanRateScatter.pdf', SubjectID, Date,NeuralInputID{1},NeuralInputID{2})),'-dpdf','-fillpage')
end
