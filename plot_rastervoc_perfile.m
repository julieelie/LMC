function plot_rastervoc_perfile(InputDataFile, OutputPath,Delay, PlotDyn, DurOrd)
%%
if nargin<3
    Delay=[3000 200];
end
if nargin<4
    PlotDyn = 1;
end

if nargin<5
    DurOrd = 0; % set to 1 to order neural responses by decreasing call duration
end
MinNumCall =8; % Minimum number of events (vocalizations) to calculate a PSTH
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
FileNameBase = sprintf('%s_%s_SS%s_%s-%s', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2});
FullDataSetFile = fullfile(OutputPath, sprintf('%s.mat', FileNameBase));
Data=load(FullDataSetFile);

% Indices of vocalization renditions produced by subject ('self') or not
IndVocP = find(contains(Data.What, 'Voc').*contains(Data.Who, 'self'));
IndVocH = find(contains(Data.What, 'Voc').*~contains(Data.Who, 'self'));
IndVocD = find(((Data.DelayBefore>=Delay(1)) + contains(Data.VocRank, 'first')).*((Data.DelayAfter>=Delay(2))+ contains(Data.VocRank, 'end')));
IndVocO = find(contains(Data.ExpType, 'O'));
IndVocF = find(contains(Data.What, 'Voc').*contains(Data.ExpType, 'F'));
IndVocPD = intersect(IndVocP, IndVocD);
IndVocPDO = intersect(IndVocPD, IndVocO);
IndVocPDF = intersect(IndVocPD, IndVocF);
IndVocHD = intersect(IndVocH, IndVocD);
IndVocHDO = intersect(IndVocHD, IndVocO);
IndVocHDF = intersect(IndVocHD, IndVocF);
IndVocFD = intersect(IndVocD,IndVocF);


%% Figure only onset with spectrograms and amplitude
% 
% Fig6bis = figure();
% % ColorLegend.name = {'Ba' 'Tr'};
% ColorLegend.name = {'Vocalization'};
% % ColorLegend.color = {[0/255 191/255 255/255]; [1 0.7 0.7]};
% ColorLegend.color = {[0.6350, 0.0780, 0.1840, 0.5]};
% % Color = ColorLegend.color{1}.*contains(Data.What, ColorLegend.name{1}) + ColorLegend.color{2}.*contains(Data.What, ColorLegend.name{2});
% Color = [0.6350, 0.0780, 0.1840, 0.5];
% if isfield(Data.KDE_onset, 'SelfVocOp')
%     ColKDE = [0.25 0.25 0.25];
%     timerasterkdeOnSpectroAmp(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPDO,Data.BioSound, Color,ColorLegend,Data.KDE_onset.SelfVocOp,ColKDE)
% end

%% Figure only onset of both hearing and producing during Free session with amplitude envelope
ColorLegend.name = {'Vocalizing', 'Hearing'};
ColorLegend.color = {[0.4940, 0.1840, 0.5560, 0.5]; [0.3010, 0.7450, 0.9330, 0.5] };
Color = ColorLegend.color{1}.*contains(Data.Who, 'self') + ColorLegend.color{2}.*(~contains(Data.Who, 'self'));
% Color = [0.6350, 0.0780, 0.1840, 0.5];
if isfield(Data.KDE_onset, 'SelfVocFr') && isfield(Data.KDE_onset, 'OthersVocFr')
    Fig10 = figure(10);
    ColKDE = [0.4940, 0.1840, 0.5560; 0.3010, 0.7450, 0.9330] ;
    timerasterkdeOnAmp2(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,[IndVocPDF; IndVocHDF],Data.BioSound, Color,ColorLegend,IndVocPDF,IndVocHDF, ColKDE )
    print(Fig10,fullfile(OutputPath,sprintf('%s_RasterVocFreeSelfOthers_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage');
    close(Fig10)
end

% %% Time Raster plot alligned to vocalization production onset/offset self vocalizations Operant + Free First voc of sequence only
% if ~isempty(IndVocPD) && ~isempty(IndVocPDO) && ~isempty(IndVocPDF) && length(IndVocPD)>MinNumCall
%     Fig1 = figure();% TrCol = [0.9290, 0.6940, 0.1250];BaCol = [1, 0, 0];
%     ColorLegend.name = {'Ba' 'Tr'};
%     ColorLegend.color = {[0/255 191/255 255/255]; [1 0.7 0.7]};
%     Color = ColorLegend.color{1}.*contains(Data.What, ColorLegend.name{1}) + ColorLegend.color{2}.*contains(Data.What, ColorLegend.name{2});
%     if isfield(Data.KDE_onset,'SelfVocAll')
%         ColKDE = [186/255 85/255 211/255];
%         timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPD,Color,ColorLegend,Data.KDE_onset.SelfVocAll,Data.KDE_offset.SelfVocAll,ColKDE,Data.RewardTime,DurOrd);
%     else
%         timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPD,Color, ColorLegend);
%     end
%     suplabel(sprintf('CALLS FROM SUBJECT O and F   %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
%     print(Fig1,fullfile(OutputPath,sprintf('%s_RasterVocSelf_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage');
% end
% 
% %% Time Raster plot alligned to vocalization perception onset/offset others vocalizations Operant + Free First voc of sequence only
% if ~isempty(IndVocHD) && ~isempty(IndVocHDO) && ~isempty(IndVocHDF) && length(IndVocHD)>MinNumCall
%     Fig2 = figure();
%     ColorLegend.name = {'Ba' 'Tr'};
%     ColorLegend.color = {[0/255 191/255 255/255]; [1 0.7 0.7]};
%     Color = ColorLegend.color{1}.*contains(Data.What, ColorLegend.name{1}) + ColorLegend.color{2}.*contains(Data.What, ColorLegend.name{2});
%     
%     if isfield(Data.KDE_onset, 'OthersVocAll')
%         ColKDE =[0.4 0.45 1];
%         timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHD,Color,ColorLegend, Data.KDE_onset.OthersVocAll,Data.KDE_offset.OthersVocAll,ColKDE,Data.RewardTime,DurOrd)
%     else
%         timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHD,Color, ColorLegend)
%     end
%     suplabel(sprintf('CALLS FROM OTHERS O and F    %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
%     print(Fig2,fullfile(OutputPath,sprintf('%s_RasterVocOthers_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
% end
% 
% %% Time Raster plot alligned to vocalization production onset/offset during Operant conditioning First voc of sequence only
% if ~isempty(IndVocPDO) && length(IndVocPDO)>MinNumCall
%     Fig6 = figure();
%     ColorLegend.name = {'Ba' 'Tr'};
%     ColorLegend.color = {[0/255 191/255 255/255]; [1 0.7 0.7]};
%     Color = ColorLegend.color{1}.*contains(Data.What, ColorLegend.name{1}) + ColorLegend.color{2}.*contains(Data.What, ColorLegend.name{2});
%     if isfield(Data.KDE_onset, 'SelfVocOp')
%         ColKDE = [186/255 85/255 211/255];
%         timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPDO,Color,ColorLegend,Data.KDE_onset.SelfVocOp,Data.KDE_offset.SelfVocOp,ColKDE,Data.RewardTime,DurOrd)
%     else
%         timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPDO,Color,ColorLegend)
%     end
%     suplabel(sprintf('CALLS FROM SUBJECT OPERANT   %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
%     print(Fig6,fullfile(OutputPath,sprintf('%s_RasterVocSelfOp_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
% end
% 
% 
% %% Time Raster plot alligned to vocalization perception onset/offset during operant conditioning First voc of sequence only
% if ~isempty(IndVocHDO) && length(IndVocHDO)>MinNumCall
%     Fig7 = figure();
%     ColorLegend.name = {'Ba' 'Tr'};
%     ColorLegend.color = {[0/255 191/255 255/255]; [1 0.7 0.7]};
%     Color = ColorLegend.color{1}.*contains(Data.What, ColorLegend.name{1}) + ColorLegend.color{2}.*contains(Data.What, ColorLegend.name{2});
%     if isfield(Data.KDE_onset,'OthersVocOp')
%         ColKDE = [0.4 0.45 1];
%         timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHDO,Color,ColorLegend,Data.KDE_onset.OthersVocOp,Data.KDE_offset.OthersVocOp,ColKDE,Data.RewardTime,DurOrd)
%     else
%         timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHDO,Color, ColorLegend)
%     end
%     suplabel(sprintf('CALLS FROM OTHERS OPERANT    %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
%     print(Fig7,fullfile(OutputPath,sprintf('%s_RasterVocOthersOp_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
% end
% 
% %% Time Raster plot alligned to vocalization production onset/offset during Free session First voc of sequence only
% if ~isempty(IndVocPDF) && length(IndVocPDF)>MinNumCall
%     Fig8 = figure();
%     ColorLegend.name = {'Ba' 'Tr'};
%     ColorLegend.color = {[0/255 191/255 255/255]; [1 0.7 0.7]};
%     Color = ColorLegend.color{1}.*contains(Data.What, ColorLegend.name{1}) + ColorLegend.color{2}.*contains(Data.What, ColorLegend.name{2});
%     if isfield(Data.KDE_onset,'SelfVocFr')
%         ColKDE = [186/255 85/255 211/255];
%         timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPDF,Color,ColorLegend,Data.KDE_onset.SelfVocFr,Data.KDE_offset.SelfVocFr,ColKDE,Data.RewardTime,DurOrd)
%     else
%         timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPDF,Color, ColorLegend)
%     end
%     suplabel(sprintf('CALLS FROM SUBJECT FREE SESSION    %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
%     print(Fig8,fullfile(OutputPath,sprintf('%s_RasterVocSelfFr_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
% end
% 
% %% Time Raster plot alligned to vocalization perception onset/offset during Free session First voc of sequence only
% if ~isempty(IndVocHDF) && length(IndVocHDF)>MinNumCall
%     Fig9 = figure();
%     ColorLegend.name = {'Ba' 'Tr'};
%     ColorLegend.color = {[0/255 191/255 255/255]; [1 0.7 0.7]};
%     Color = ColorLegend.color{1}.*contains(Data.What, ColorLegend.name{1}) + ColorLegend.color{2}.*contains(Data.What, ColorLegend.name{2});
%     if isfield(Data.KDE_onset,'OthersVocFr')
%         ColKDE = [0.4 0.45 1];
%         timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHDF,Color,ColorLegend, Data.KDE_onset.OthersVocFr,Data.KDE_offset.OthersVocFr,ColKDE,Data.RewardTime,DurOrd)
%     else
%         timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHDF,Color, ColorLegend)
%     end
%     suplabel(sprintf('CALLS FROM OTHERS FREE SESSION    %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
%     print(Fig9,fullfile(OutputPath,sprintf('%s_RasterVocOthersFr_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
% end

%% DYNAMIC RASTERS
if PlotDyn
    %% Dynamic raster plot All Calls
    % count the number of spike during vocalization, if no spike no plot
    Ind_SAT = nan(length(IndVocP),1);
    for vv=1:length(IndVocP)
        SAT = Data.SpikesArrivalTimes_Behav{IndVocP(vv)};
        Ind_SAT(vv) = sum((SAT>0).*(SAT<=Data.Duration(IndVocP(vv))));
    end
    if any(sum(Ind_SAT))
        fprintf(1,'All Calls raster plot\n')
    %     for vv=1:length(IndVocP)
    %         fprintf(1,'All calls ratser: Voc %d/%d\n',vv, length(IndVocP))
    %         Fig3=figure(10);
    %         if vv==1
    %             Legend=1;
    %         else
    %             Legend=0;
    %         end
    %         hold on
    %         plotCallDynamic(Data.BioSound{IndVocP(vv),1}, Data.BioSound{IndVocP(vv),2},[],Legend)
    %     end
    % Just plot one call for the legend
    plotCallDynamic(Data.BioSound{IndVocP(1),1}, Data.BioSound{IndVocP(1),2},[],1)
        for vv=1:length(IndVocP)
            fprintf(1,'All calls ratser: Voc %d/%d\n',vv, length(IndVocP))
            Fig3=figure(10);
            SAT = Data.SpikesArrivalTimes_Behav{IndVocP(vv)};
            if any(Ind_SAT(vv))
                plotCallDynamic(Data.BioSound{IndVocP(vv),1}, Data.BioSound{IndVocP(vv),2},SAT(Ind_SAT(vv)),0);
            end
        end
        suplabel(sprintf('CALLS FROM SUBJECT O and F   %s on %s Raster Tetrode %s S%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
        print(Fig3,fullfile(OutputPath,sprintf('%s_DynRasterVocProd_all.pdf', FileNameBase)),'-dpdf','-fillpage')
    end
    %% Dynamic raster plot Trills
    IndVocPTr = find(contains(Data.What, 'Tr') .*contains(Data.Who, 'self'));
    if ~isempty(IndVocPTr)
        Ind_SAT = nan(length(IndVocPTr),1);
        for vv=1:length(IndVocPTr)
            SAT = Data.SpikesArrivalTimes_Behav{IndVocPTr(vv)};
            Ind_SAT(vv) = sum((SAT>0).*(SAT<=Data.Duration(IndVocPTr(vv))));
        end
        if any(sum(Ind_SAT))
            fprintf(1,'Trill raster plot\n')
            for vv=1:length(IndVocPTr)
                fprintf(1,'Trill raster Voc %d/%d\n',vv, length(IndVocPTr))
                Fig4=figure(11);
                if vv==1
                    Legend=1;
                else
                    Legend=0;
                end
                hold on
                plotCallDynamic(Data.BioSound{IndVocPTr(vv),1}, Data.BioSound{IndVocPTr(vv),2},[],Legend)
            end
            for vv=1:length(IndVocPTr)
                fprintf(1,'Trill raster Neural data Voc %d/%d\n',vv, length(IndVocPTr))
                Fig4=figure(11);
                SAT = Data.SpikesArrivalTimes_Behav{IndVocPTr(vv)};
                if any(Ind_SAT(vv))
                    plotCallDynamic(Data.BioSound{IndVocPTr(vv),1}, Data.BioSound{IndVocPTr(vv),2},SAT(Ind_SAT(vv)),0);
                end
            end
            suplabel(sprintf('TRILLS FROM SUBJECT O and F   %s on %s Raster Tetrode %s S%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
            print(Fig4,fullfile(OutputPath,sprintf('%s_DynRasterVocProd_Tr.pdf', FileNameBase)),'-dpdf','-fillpage')
        end
    end
    %% Dynamic raster plot Barks
    IndVocPBa = find(contains(Data.What, 'Ba').*contains(Data.Who, 'self'));
    if ~isempty(IndVocPBa)
        Ind_SAT = nan(length(IndVocPBa),1);
        for vv=1:length(IndVocPBa)
            SAT = Data.SpikesArrivalTimes_Behav{IndVocPBa(vv)};
            Ind_SAT(vv) = sum((SAT>0).*(SAT<=Data.Duration(IndVocPBa(vv))));
        end
        if any(sum(Ind_SAT))
            fprintf(1,'Bark raster plot\n')
            for vv=1:length(IndVocPBa)
                fprintf(1,'Bark raster Voc %d/%d\n',vv, length(IndVocPBa))
                Fig5=figure(12);
                if vv==1
                    Legend=1;
                else
                    Legend=0;
                end
                hold on
                plotCallDynamic(Data.BioSound{IndVocPBa(vv),1}, Data.BioSound{IndVocPBa(vv),2},[],Legend)
            end
            for vv=1:length(IndVocPBa)
                fprintf(1,'Bark raster neural data Voc %d/%d\n',vv, length(IndVocPBa))
                Fig5=figure(12);
                SAT = Data.SpikesArrivalTimes_Behav{IndVocPBa(vv)};
                if any(Ind_SAT(vv))
                    plotCallDynamic(Data.BioSound{IndVocPBa(vv),1}, Data.BioSound{IndVocPBa(vv),2},SAT(Ind_SAT(vv)),0);
                end
            end
            suplabel(sprintf('BARKS FROM SUBJECT O and F   %s on %s Raster Tetrode %s S%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
            print(Fig5,fullfile(OutputPath,sprintf('%s_DynRasterVocProd_Ba.pdf', FileNameBase)),'-dpdf','-fillpage')
        end
    end
end
close all

%% INTERNAL FUNCTION
    function timeraster(SpikesArrivalTimes,Duration,Delay,Indices, Color, ColorLegend)
        % We want to plot data with increasing duration of
        % vocalizations.
        [~,IDur] = sort(Duration(Indices));

        ss1=subplot(1,2,1);

        % First alligned to vocalization onset
        for oo=1:length(Indices)
            cc=IDur(oo);
            hold on
            plot([0 Duration(Indices(cc))], oo-[0.5 0.5], '-','LineWidth',250/length(Indices),'Color', Color(Indices(cc),:)) % vocalization
            Sat = SpikesArrivalTimes{Indices(cc)};
            SpikeInd = find((Sat>-Delay(1)).*(Sat<(Duration(Indices(cc))+Delay(2))));
            if ~isempty(SpikeInd)
                Sat = Sat(SpikeInd);
                for spike=1:length(Sat)
                    hold on
                    plot(Sat(spike)*ones(2,1), oo-[0.9 0.1], 'k-', 'LineWidth',2)
                end
            end
            hold on
        end
        XLIM = [-Delay(1) max(Duration(Indices))+Delay(2)];
        xlabel('Time centered at vocalization onset (ms)')
        ylim([0 length(Indices)+1])
        xlim(XLIM)
        ylabel('Vocalization renditions')
        title(ss1, ['\fontsize{16} {\color[rgb]' sprintf('{%.1f %.1f %.1f}%s',ColorLegend.color{1}, ColorLegend.name{1}) ' \color[rgb]' sprintf('{%.1f %.1f %.1f}%s}', ColorLegend.color{2}, ColorLegend.name{2})]);
        hold off
        
        
        % then alligned to vocalization offset
        ss2= subplot(1,2,2);

        for oo=1:length(Indices)
            cc=IDur(oo);
            hold on
            plot([-Duration(Indices(cc)) 0], oo-[0.5 0.5], '-','LineWidth',250/length(Indices),'Color', Color(Indices(cc),:)) % vocalization
            Sat = SpikesArrivalTimes{Indices(cc)};
            SpikeInd = find((Sat>-Delay(1)).*(Sat<(Duration(Indices(cc))+Delay(2))));
            if ~isempty(SpikeInd)
                Sat = Sat(SpikeInd)- Duration(Indices(cc));
                for spike=1:length(Sat)
                    hold on
                    plot(Sat(spike)*ones(2,1), oo-[0.9 0.1], 'k-', 'LineWidth',2)
                end
            end
            hold on
        end
        XLIM = [-Delay(1)-max(Duration(Indices)) Delay(2)];
        xlabel('Time centered at vocalization offset (ms)')
        ylim([0 length(Indices)+1])
        xlim(XLIM)
        ylabel('Vocalization renditions')
        title(ss2, ['\fontsize{16} {\color[rgb]' sprintf('{%.1f %.1f %.1f}%s',ColorLegend.color{1}, ColorLegend.name{1}) ' \color[rgb]' sprintf('{%.1f %.1f %.1f}%s}', ColorLegend.color{2}, ColorLegend.name{2})]);
        hold off
    end


function timerasterkde(SpikesArrivalTimes,Duration,Delay,Indices, Color, ColorLegend, Dat1, Dat2, Col1, RewardTime,DurOrd)
    if nargin<9
        RewardTime = nan(length(Duration),1);
    end
    if nargin<10
        DurOrd=0;
    end
    if DurOrd
        % We want to plot data with increasing duration of
        % vocalizations.
        [~,IDur] = sort(Duration(Indices));
    else
        IDur = 1:length(Indices);
    end
    
        ss1=subplot(4,2,[1 3 5]);

        % First alligned to vocalization onset
        for oo=1:length(Indices)
            cc=IDur(oo);
            hold on
            % plot vocalization timing
            plot([0 Duration(Indices(cc))], oo-[0.5 0.5], '-','LineWidth',250/length(Indices),'Color', Color(Indices(cc),:)) % vocalization
            % plot reward time
            if ~isinf(RewardTime(Indices(cc))) && ~isnan(RewardTime(Indices(cc)))
                plot(RewardTime(Indices(cc)), oo-0.5,'o','MarkerSize',6, 'MarkerFaceColor',[1 0.85 0.275],'MarkerEdgeColor',[1 0.85 0.275])
                hold on
            end
            % plot spikes
            Sat = SpikesArrivalTimes{Indices(cc)};
            if Delay(2)<5000
                SpikeInd = find((Sat>-Delay(1)).*(Sat<(max(Duration(Indices))+Delay(2))));
            else
                SpikeInd = find((Sat>-Delay(1)).*(Sat<(Duration(Indices(cc))+Delay(2))));
            end
            if ~isempty(SpikeInd)
                Sat = Sat(SpikeInd);
                for spike=1:length(Sat)
                    hold on
                    plot(Sat(spike)*ones(2,1), oo-[0.9 0.1], 'k-', 'LineWidth',1.5)
                end
            end
            hold on
            
        end
        XLIM = [-Delay(1) max(Duration(Indices))+Delay(2)];
        xlabel('Time centered at vocalization onset (ms)')
        ylim([0 length(Indices)+1])
        xlim(XLIM)
        ylabel('Vocalization renditions')
        title(ss1, ['\fontsize{16} {\color[rgb]' sprintf('{%.1f %.1f %.1f}%s',ColorLegend.color{1}, ColorLegend.name{1}) ' \color[rgb]' sprintf('{%.1f %.1f %.1f}%s}', ColorLegend.color{2}, ColorLegend.name{2})]);
        hold off
        
        subplot(4,2,7)
        plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
        hold on
        shadedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim(XLIM)
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        
        % then alligned to vocalization offset
        ss2=subplot(4,2,[2 4 6]);

        for oo=1:length(Indices)
            cc=IDur(oo);
            hold on
            % plot call duration
            plot([-Duration(Indices(cc)) 0], oo-[0.5 0.5], '-','LineWidth',250/length(Indices),'Color', Color(Indices(cc),:)) % vocalization
            % plot reward time
            if ~isinf(RewardTime(Indices(cc))) && ~isnan(RewardTime(Indices(cc)))
                plot(RewardTime(Indices(cc))- Duration(Indices(cc)), oo-0.5,'o','MarkerSize',8, 'MarkerFaceColor',[1 0.85 0.275],'MarkerEdgeColor',[1 0.85 0.275])
                hold on
            end
            % plot spikes
            Sat = SpikesArrivalTimes{Indices(cc)};
            if Delay(1)<5000
                SpikeInd = find((Sat>(-Delay(1)-max(Duration(Indices))+Duration(Indices(cc)))).*(Sat<(Duration(Indices(cc))+Delay(2))));
            else
                SpikeInd = find((Sat>-Delay(1)).*(Sat<(Duration(Indices(cc))+Delay(2))));
            end
            if ~isempty(SpikeInd)
                Sat = Sat(SpikeInd)- Duration(Indices(cc));
                for spike=1:length(Sat)
                    hold on
                    plot(Sat(spike)*ones(2,1), oo-[0.9 0.1], 'k-', 'LineWidth',2)
                end
            end
            hold on
            
        end
        XLIM = [-Delay(1)-max(Duration(Indices)) Delay(2)];
        xlabel('Time centered at vocalization offset (ms)')
        ylim([0 length(Indices)+1])
        xlim(XLIM)
        ylabel('Vocalization renditions')
        title(ss2, ['\fontsize{16} {\color[rgb]' sprintf('{%.1f %.1f %.1f}%s',ColorLegend.color{1}, ColorLegend.name{1}) ' \color[rgb]' sprintf('{%.1f %.1f %.1f}%s}', ColorLegend.color{2}, ColorLegend.name{2})]);
        hold off
        
        subplot(4,2,8)
        plot(Dat2(2,:),Dat2(1,:),'-','Color',Col1,'LineWidth',2)
        hold on
        shadedErrorBar(Dat2(2,:),Dat2(1,:),Dat2(3:4,:),{'-','Color',Col1, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim(XLIM)
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
    end


    function plotCallDynamic(BiosoundRaw, BiosoundPiezo,SAT,Legend)
        SpikeMarkerSize = 5;
        MAP = colormap('hot');
        MAP = flip(MAP(1:10:128,:));
        Span = 9;% Span is an unevennumber. smooth has a default span of 5 points = 5ms However end points are unchanged...
        HalfSpan = (Span-1)/2;
        % Plot the pitch saliency vs amplitude on microphone
        subplot(4,1,1)
        Saliency = mysmooth(double(BiosoundRaw.sal), Span);
        TimeSound = double(BiosoundRaw.to)*1000;
        TimeSound_norm = TimeSound./max(TimeSound);
        cmap = colormap('jet');
        ncolors = length(cmap);
        nx = length(Saliency);
        if isempty(SAT)
            for ii=HalfSpan:nx-HalfSpan
                segcolor = cmap(fix((TimeSound_norm(ii)+TimeSound_norm(ii+1))*ncolors./3)+1,:);
                plot([Saliency(ii), Saliency(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
                hold on;
            end
        else
            for ii=HalfSpan:nx-HalfSpan
                Nspike = sum((SAT>=TimeSound(ii)) .* (SAT<TimeSound(ii+1)));
                if Nspike
                    SpikeMarkerColor = MAP(Nspike,:);
                    plot([Saliency(ii), Saliency(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)],'o','MarkerSize',SpikeMarkerSize,'MarkerFaceColor', SpikeMarkerColor, 'MarkerEdgeColor',SpikeMarkerColor);
                end
                hold on
            end
        end
        set(gca,'XLim',[0 1]);
        if Legend
            xlabel('Pitch Saliency')
            ylabel('Amplitude')
        end
        
        % Plot the difference of formants (Mic data) vs sound amplitude (Mic
        % Data)
        subplot(4,1,2)
        SoundSpeed = 350;
        F1 = double(BiosoundRaw.F1);
        F2 = double(BiosoundRaw.F2);
        FormantDisp = mysmooth(SoundSpeed./(2*(F2 - F1))*1000, Span);
        nx = length(FormantDisp);
        
        if isempty(SAT)
            for ii=HalfSpan:nx-HalfSpan
                segcolor = cmap(fix((TimeSound_norm(ii)+TimeSound_norm(ii+1))*ncolors./3)+1,:);
                plot([FormantDisp(ii), FormantDisp(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
                hold on;
            end
        else
            for ii=HalfSpan:nx-HalfSpan
                Nspike = sum((SAT>=TimeSound(ii)) .* (SAT<TimeSound(ii+1)));
                if Nspike
                    SpikeMarkerColor = MAP(Nspike,:);
                    plot([FormantDisp(ii), FormantDisp(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)],'o','MarkerSize',SpikeMarkerSize,'MarkerFaceColor', SpikeMarkerColor, 'MarkerEdgeColor',SpikeMarkerColor);
                end
                hold on
            end
        end
        set(gca,'XLim',[10 150])
        if Legend
            xlabel('1/Formant disp (vocal tract length (mm))')
            ylabel('Amplitude')
        end
        
        % Plot the amplitude (Mic data) vs fundamental (Piezo
        % Data)
        subplot(4,1,3)
        SoundFund = mysmooth(double(BiosoundPiezo.f0), Span);
        if ~isempty(SoundFund)
            if isempty(SAT)
                for ii=HalfSpan:nx-HalfSpan
                    segcolor = cmap(fix((TimeSound_norm(ii)+TimeSound_norm(ii+1))*ncolors./3)+1,:);
                    plot([SoundFund(ii), SoundFund(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
                    hold on;
                end
            else
                for ii=HalfSpan:nx-HalfSpan
                    Nspike = sum((SAT>=TimeSound(ii)) .* (SAT<TimeSound(ii+1)));
                    if Nspike
                        SpikeMarkerColor = MAP(Nspike,:);
                        plot([SoundFund(ii), SoundFund(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)],'o','MarkerSize',SpikeMarkerSize,'MarkerFaceColor', SpikeMarkerColor, 'MarkerEdgeColor',SpikeMarkerColor);
                    end
                    hold on
                end
            end
            if Legend
                ylabel('Amplitude')
                xlabel('Fundamental (Hz)')
            end
            set(gca,'XLim',[200 3000])
        end
        
        % Plot the Amplitude (Mic data) vs SpectralMean (Mic
        % Data)
        subplot(4,1,4)
        SoundSpecMean = mysmooth(double(BiosoundRaw.SpectralMean), Span);
        if ~isempty(SoundSpecMean)
            if isempty(SAT)
                for ii=HalfSpan:nx-HalfSpan
                    segcolor = cmap(fix((TimeSound_norm(ii)+TimeSound_norm(ii+1))*ncolors./3)+1,:);
                    plot([SoundSpecMean(ii), SoundSpecMean(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
                    hold on;
                end
            else
                for ii=HalfSpan:nx-HalfSpan
                    Nspike = sum((SAT>=TimeSound(ii)) .* (SAT<TimeSound(ii+1)));
                    if Nspike
                        SpikeMarkerColor = MAP(Nspike,:);
                        plot([SoundSpecMean(ii), SoundSpecMean(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)],'o','MarkerSize',SpikeMarkerSize,'MarkerFaceColor', SpikeMarkerColor, 'MarkerEdgeColor',SpikeMarkerColor);
                    end
                    hold on
                end
            end
            if Legend
                ylabel('Amplitude')
                xlabel('Spectral Mean (Hz)')
            end
            set(gca,'XLim',[25000 30000])
        end
        
        %         % Plot the Amplitude (Mic data) vs Spectral Max (Piezo
        %         % Data)
        %         subplot(5,1,5)
        %         SoundSpecMax = mysmooth(double(BiosoundPiezo.SpectralMax), Span);
        %         if ~isempty(SoundSpecMax)
        %             for ii=HalfSpan:nx-HalfSpan
        %                 segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
        %                 plot([SoundSpecMax(ii), SoundSpecMax(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
        %                 hold on;
        %             end
        %             ylabel('Amplitude')
        %             xlabel(sprintf('Spectral Max (Hz), %.1f Hz', nanmean(double(BiosoundPiezo.SpectralMax))))
        %             set(gca,'XLim',[0 10000])
        %         end
        
    end


    function outyy = mysmooth(yy,Span)
        if nargin<2
            Span = 5;
        end
        outyy=nan(size(yy));
        for ii=1:length(yy)
            if ii==1 || ii==length(yy)
                outyy(ii) = yy(ii);
            elseif ii<((Span-1)/2)
                HalfSpan = ii-1;
                outyy(ii) = nanmean(yy(1:(ii+HalfSpan)));
            elseif (length(yy)-ii) < ((Span-1)/2)
                HalfSpan = length(yy)-ii;
                outyy(ii) = nanmean(yy((ii-HalfSpan):end));
            else
                outyy(ii) = nanmean(yy((ii-HalfSpan):(ii+HalfSpan)));
            end
        end
    end


function timerasterkdeOnSpectroAmp(SpikesArrivalTimes,Duration,Delay,Indices, BioSound, Color, ColorLegend, Dat1, Col1, RewardTime,DurOrd)
    if nargin<10
        RewardTime = nan(length(Duration),1);
    end
    if nargin<11
        DurOrd=0;
    end
    if DurOrd
        % We want to plot data with increasing duration of
        % vocalizations.
        [~,IDur] = sort(Duration(Indices));
    else
        IDur = 1:length(Indices);
    end
    
    % Plot spectrogram of the longest call from Microphone
    ss1=subplot(7,1,1);
    DBNOISE =60;
    f_low = 0;
    LongestVoc = find(Duration==max(Duration(Indices)));
    LongestVoc = LongestVoc(1);
    yyaxis left
    logB = BioSound{LongestVoc,1}.spectro;
    maxB = max(max(logB));
    minB = maxB-DBNOISE;
    imagesc(double(BioSound{LongestVoc,1}.to)*1000,double(BioSound{LongestVoc,1}.fo),logB);          % to is in seconds
    axis xy;
    caxis('manual');
    caxis([minB maxB]);
    cmap = spec_cmap();
    colormap(cmap);
    %         colorbar()
    v_axis = axis;
    v_axis(3)=f_low;
    v_axis(4)=50000;
    axis(v_axis);
    xlabel('time (ms)'), ylabel('Frequency');
    XLIM = [-Delay(1) max(Duration(Indices))+Delay(2)];
    xlim(XLIM)
    
    % Plot spectrogram of the longest call from Piezo
    ss2=subplot(7,1,2);
    DBNOISE =60;
    f_low = 0;
    yyaxis left
    logB = BioSound{LongestVoc,2}.spectro;
    maxB = max(max(logB));
    minB = maxB-DBNOISE;
    imagesc(double(BioSound{LongestVoc,2}.to)*1000,double(BioSound{LongestVoc,2}.fo),logB);          % to is in seconds
    axis xy;
    caxis('manual');
    caxis([minB maxB]);
    cmap = spec_cmap();
    colormap(cmap);
    %         colorbar()
    v_axis = axis;
    v_axis(3)=f_low;
    v_axis(4)=10000;
    axis(v_axis);
    xlabel('time (ms)'), ylabel('Frequency');
    XLIM = [-Delay(1) max(Duration(Indices))+Delay(2)];
    xlim(XLIM)        
    
    % Plot spikes and vocalization spots alligned to vocalization onset
    % and gather amplitude data
    XAmp_Mat = zeros(length(Indices), Duration(LongestVoc));
    YSpike_Mat = zeros(length(Indices), Duration(LongestVoc) + sum(Delay));
    ss3=subplot(7,1,[3 4 5]);
    for oo=1:length(Indices)
        cc=IDur(oo);
        XAmp_Mat(oo,round(BioSound{Indices(oo),2}.tAmp .*1000)+1) = BioSound{Indices(oo),2}.amp;
        
        hold on
        % plot vocalization timing
        if size(Color,1)==1
            plot([0 Duration(Indices(cc))], oo-[0.5 0.5], '-','LineWidth',250/length(Indices),'Color', Color)
        else
            plot([0 Duration(Indices(cc))], oo-[0.5 0.5], '-','LineWidth',250/length(Indices),'Color', Color(Indices(cc),:)) % vocalization
        end
        
        % plot reward time
        if ~isinf(RewardTime(Indices(cc))) && ~isnan(RewardTime(Indices(cc)))
            plot(RewardTime(Indices(cc)), oo-0.5,'o','MarkerSize',6, 'MarkerFaceColor',[1 0.85 0.275],'MarkerEdgeColor',[1 0.85 0.275])
            hold on
        end
        % plot spikes
        Sat = SpikesArrivalTimes{Indices(cc)};
        if Delay(2)<5000
            SpikeInd = find((Sat>-Delay(1)).*(Sat<(max(Duration(Indices))+Delay(2))));
        else
            SpikeInd = find((Sat>-Delay(1)).*(Sat<(Duration(Indices(cc))+Delay(2))));
        end
        if ~isempty(SpikeInd)
            Sat = Sat(SpikeInd);
            for spike=1:length(Sat)
                hold on
                plot(Sat(spike)*ones(2,1), oo-[0.9 0.1], 'k-', 'LineWidth',1.5)
                YSpike_Mat(oo, ceil(Sat(spike)+Delay(1)))=1;
            end
        end
        hold on
        
    end
    XLIM = [-Delay(1) max(Duration(Indices))+Delay(2)];
    xlabel('Time centered at vocalization onset (ms)')
    ylim([0 length(Indices)+1])
    xlim(XLIM)
    ylabel('Vocalization renditions')
    if length(ColorLegend.color)>1
        title(ss1, ['\fontsize{16} {\color[rgb]' sprintf('{%.1f %.1f %.1f}%s',ColorLegend.color{1}, ColorLegend.name{1}) ' \color[rgb]' sprintf('{%.1f %.1f %.1f}%s}', ColorLegend.color{2}, ColorLegend.name{2})]);
    else
%         title(ss1, ['\fontsize{16} {\color[rgb]' sprintf('{%.1f %.1f %.1f}%s',ColorLegend.color{1}(1:3), ColorLegend.name{1})]);
    end
    hold off
    
    % Plot the KDE and on top the average amplitude of vocalizations
    TR=2;
    Overlap = 0;
    [YPerStim, YPerStimt,FS] = get_y_4Coherence(SpikesArrivalTimes(Indices), Duration(Indices),Delay,TR,Overlap);
    
    subplot(7,1,[6 7])
    yyaxis left
    cla
    AvAmp = [zeros(1, Delay(1)) mean(XAmp_Mat,1) zeros(1, Delay(2))];
    Ampx = -(Delay(1)):1:max(Duration(Indices))+Delay(2)-1;
    StdAmp = [zeros(1, Delay(1)) std(XAmp_Mat,0,1)./(size(XAmp_Mat,1))^0.5 zeros(1, Delay(2))];
    shadedErrorBar(Ampx,AvAmp,StdAmp,{'-','Color',[0.6350, 0.0780, 0.1840, 0.7], 'LineWidth',2})
    ylabel('Vocalization Amplitude')
    
    yyaxis right
    cla
    shadedErrorBar(YPerStimt{1}, mean(cell2mat(YPerStim'),1).*1000, 1000*std(cell2mat(YPerStim'),0,1)./(length(YPerStim)).^0.5, {'-','Color',Col1, 'LineWidth',2})
%     plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
%     hold on
%     shadedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
    hold on
    VL = vline(0,':k');
    VL.LineWidth = 2;
    hold off
    xlim(XLIM)
    xlabel('Time (ms)')
    ylabel('Rate (Hz)')
    
    
    
    figure()
    [z,p,k] = butter(6,60/(FS/2),'low');
    sos_lowY = zp2sos(z,p,k);
    [z,p,k] = butter(6,60/(1000/2),'low');
    sos_lowX = zp2sos(z,p,k);
    YPerStimLow = YPerStim;
    XAmp_MatLow = XAmp_Mat;
    for Stim=1:length(YPerStim)
        YPerStimLow{Stim} = filtfilt(sos_lowY,1,YPerStim{Stim});
        XAmp_MatLow(Stim,:) = filtfilt(sos_lowX,1,XAmp_Mat(Stim,:));
    end
    yyaxis left
    cla
    AvAmp = [zeros(1, Delay(1)) mean(XAmp_MatLow,1) zeros(1, Delay(2))];
    Ampx = -(Delay(1)):1:max(Duration(Indices))+Delay(2)-1;
    StdAmp = [zeros(1, Delay(1)) std(XAmp_MatLow,0,1)./(size(XAmp_MatLow,1))^0.5 zeros(1, Delay(2))];
    shadedErrorBar(Ampx,AvAmp,StdAmp,{'-','Color',[0.6350, 0.0780, 0.1840, 0.7], 'LineWidth',2})
    ylabel('Vocalization Amplitude')
    
    yyaxis right
    cla
    shadedErrorBar(YPerStimt{1}, mean(cell2mat(YPerStimLow'),1).*1000, 1000*std(cell2mat(YPerStimLow'),0,1)./(length(YPerStimLow)).^0.5, {'-','Color',Col1, 'LineWidth',2})
%     plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
%     hold on
%     shadedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
    hold on
    VL = vline(0,':k');
    VL.LineWidth = 2;
    hold off
    xlim(XLIM)
    xlabel('Time (ms)')
    ylabel('Rate (Hz)')
    
    % Calculate the Amplitude MRF
    AmpWindow = [0 50];
    NStims=length(Indices);
    AllAmp = [zeros(NStims, Delay(1)) XAmp_MatLow zeros(NStims, Delay(2))];
    YPerStimLowMat = cell2mat(YPerStimLow');
    StepIndices = find((YPerStimt{1}>(-Delay(1)-AmpWindow(1))).*(YPerStimt{1}<(max(YPerStimt{1})-AmpWindow(2))));
    Nsteps = length(StepIndices);
    ConvMat = nan(Nsteps*NStims,length(AmpWindow(1):AmpWindow(2)));
    
    
    for oo=1:NStims
        for step = 1:Nsteps
            XAmp_startInd = YPerStimt{1}(StepIndices(step)) + Delay(1);
            XAmp_local = AllAmp(oo, XAmp_startInd+(AmpWindow(1):AmpWindow(2)));
            ConvMat((oo-1)*Nsteps+step,:)= XAmp_local.* YPerStimLow{oo}(StepIndices(step));
        end
    end
    figure()
    shadedErrorBar(AmpWindow(1) : AmpWindow(2),mean(ConvMat),std(ConvMat,0,1)./(size(ConvMat,1))^0.5,{'-','Color',[0.6350, 0.0780, 0.1840, 0.7], 'LineWidth',2})
end


    function timerasterkdeOnAmp2(SpikesArrivalTimes,Duration,Delay,Indices, BioSound, Color, ColorLegend, Indices1,Indices2, ColKDE, RewardTime,DurOrd)
    if nargin<11 || isnan(RewardTime)
        RewardTime = nan(length(Duration),1);
    end
    if nargin<12
        DurOrd=0;
    end
    if DurOrd==1
        % We want to plot data with increasing duration of
        % vocalizations within each category.
        [~,IDur1] = sort(Duration(Indices1));
        [~,IDur2] = sort(Duration(Indices2));
        IDur = [IDur1 IDur2+length(IDur1)];
    else
        IDur = 1:length(Indices);
    end
    
%     % Plot spectrogram of the longest call from Microphone
%     ss1=subplot(7,1,1);
%     DBNOISE =60;
%     f_low = 0;
    LongestVoc = find(Duration==max(Duration(Indices)));
    LongestVoc = LongestVoc(1);
% %     yyaxis left
%     logB = BioSound{LongestVoc,1}.spectro;
%     maxB = max(max(logB));
%     minB = maxB-DBNOISE;
%     imagesc(double(BioSound{LongestVoc,1}.to)*1000,double(BioSound{LongestVoc,1}.fo),logB);          % to is in seconds
%     axis xy;
%     caxis('manual');
%     caxis([minB maxB]);
%     cmap = spec_cmap();
%     colormap(cmap);
%     %         colorbar()
%     v_axis = axis;
%     v_axis(3)=f_low;
%     v_axis(4)=50000;
%     axis(v_axis);
%     xlabel('time (ms)'), ylabel('Frequency');
%     XLIM = [-Delay(1) max(Duration(Indices))+Delay(2)];
%     xlim(XLIM)
%     ss1.Children.Parent.YColor = ss1.Children.Parent.XColor;
%     ss1.YAxis(2).Visible = 'off';
    
    % Plot spectrogram of the longest call from Piezo 
    ss1=subplot(7,1,1);
    DBNOISE =60;
    f_low = 0;
%     yyaxis left
    logB = BioSound{LongestVoc,2}.spectro;
    maxB = max(max(logB));
    minB = maxB-DBNOISE;
    imagesc(double(BioSound{LongestVoc,2}.to)*1000,double(BioSound{LongestVoc,2}.fo),logB);          % to is in seconds
    axis xy;
    caxis('manual');
    caxis([minB maxB]);
    cmap = spec_cmap();
    colormap(cmap);
    %         colorbar()
    v_axis = axis;
    v_axis(3)=f_low;
    v_axis(4)=5000;
    axis(v_axis);
%     xlabel('time (ms)');
    ylabel('Frequency');
    XLIM = [-Delay(1) max(Duration(Indices))+Delay(2)];
    xlim(XLIM) 
    ss1.Children.Parent.YColor = ss1.Children.Parent.XColor;
    ss1.Box = 'off';
%     ss2.YAxis(2).Visible = 'off';
    
    % Plot spikes and vocalization spots alligned to vocalization onset
    % and gather amplitude data
    XAmp_Mat = zeros(length(Indices), Duration(LongestVoc));
    YSpike_Mat = zeros(length(Indices), Duration(LongestVoc) + sum(Delay));
    ss3=subplot(7,1,[4 5 6 7]);
    for oo=1:length(Indices)
        cc=IDur(oo);
        XAmp_Mat(oo,round(BioSound{Indices(cc),2}.tAmp .*1000)+1) = BioSound{Indices(cc),2}.amp;
        
        hold on
        % plot vocalization timing
        if size(Color,1)==1
            plot([0 Duration(Indices(cc))], oo-[0.5 0.5], '-','LineWidth',250/length(Indices),'Color', Color)
        else
            plot([0 Duration(Indices(cc))], oo-[0.5 0.5], '-','LineWidth',250/length(Indices),'Color', Color(Indices(cc),:)) % vocalization
%             ColSpike = Color(Indices(cc),1:3);
            ColSpike = 'k';
        end
        
        % plot reward time
        if ~isinf(RewardTime(Indices(cc))) && ~isnan(RewardTime(Indices(cc)))
            plot(RewardTime(Indices(cc)), oo-0.5,'o','MarkerSize',6, 'MarkerFaceColor',[1 0.85 0.275],'MarkerEdgeColor',[1 0.85 0.275])
            hold on
        end
        % plot spikes
        Sat = SpikesArrivalTimes{Indices(cc)};
        if Delay(2)<5000
            SpikeInd = find((Sat>-Delay(1)).*(Sat<(max(Duration(Indices))+Delay(2))));
        else
            SpikeInd = find((Sat>-Delay(1)).*(Sat<(Duration(Indices(cc))+Delay(2))));
        end
        if ~isempty(SpikeInd)
            Sat = Sat(SpikeInd);
            for spike=1:length(Sat)
                hold on
                plot(Sat(spike)*ones(2,1), oo-[0.9 0.1], 'Color', ColSpike, 'LineWidth',1.5)
                YSpike_Mat(oo, ceil(Sat(spike)+Delay(1)))=1;
            end
        end
        hold on
        
    end
    XLIM = [-Delay(1) max(Duration(Indices))+Delay(2)];
    xlabel('Time centered at vocalization onset (ms)')
    ylim([0 length(Indices)+1])
    xlim(XLIM)
    ylabel('Vocalization renditions')
    if length(ColorLegend.color)>1
        title(ss1, ['\fontsize{16} {\color[rgb]' sprintf('{%.1f %.1f %.1f}%s',ColorLegend.color{1}(1:3), ColorLegend.name{1}) ' \color[rgb]' sprintf('{%.1f %.1f %.1f}%s}', ColorLegend.color{2}(1:3), ColorLegend.name{2})]);
    else
%         title(ss1, ['\fontsize{16} {\color[rgb]' sprintf('{%.1f %.1f %.1f}%s',ColorLegend.color{1}(1:3), ColorLegend.name{1})]);
    end
    hold off
    
    
    % Plot the KDE and on top of the spectrogram the average amplitude of vocalizations
    TR=10;
    Overlap = 0;
    [YPerStim1, YPerStimt1,FS1] = get_y_4Coherence(SpikesArrivalTimes(Indices1), Duration(Indices1),Delay,TR,Overlap);
    [YPerStim2, YPerStimt2,FS2] = get_y_4Coherence(SpikesArrivalTimes(Indices2), Duration(Indices2),Delay,TR,Overlap);
    
    ss1 = subplot(7,1,1);
    yyaxis right
    cla
%     ColorAmp = [0.6350, 0.0780, 0.1840, 0.7];
    ColorAmp = [0.3, 0.3, 0.3];
    AvAmp = [zeros(1, Delay(1)) mean(XAmp_Mat,1) zeros(1, Delay(2))];
    Ampx = -(Delay(1)):1:max(Duration(Indices))+Delay(2)-1;
    StdAmp = [zeros(1, Delay(1)) std(XAmp_Mat,0,1)./(size(XAmp_Mat,1))^0.5 zeros(1, Delay(2))];
    shadedErrorBar(Ampx,AvAmp,StdAmp,{'--','Color',ColorAmp, 'LineWidth',2})
    ylabel('Average Amplitude')
    ss1.YAxis(2).Color = ColorAmp(1:3);
    
    ss4 = subplot(7,1,[2 3]);
    cla
    shadedErrorBar(YPerStimt1{1}, mean(cell2mat(YPerStim1'),1).*1000, 1000*std(cell2mat(YPerStim1'),0,1)./(length(YPerStim1)).^0.5, {'-','Color',ColKDE(1,:), 'LineWidth',2})
%     plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
%     hold on
%     shadedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
    hold on
    shadedErrorBar(YPerStimt2{1}, mean(cell2mat(YPerStim2'),1).*1000, 1000*std(cell2mat(YPerStim2'),0,1)./(length(YPerStim2)).^0.5, {'-','Color',ColKDE(2,:), 'LineWidth',2})
    hold on
    VL = vline(0,':k');
    VL.LineWidth = 2;
    hold off
    xlim(XLIM)
%     xlabel('Time (ms)')
    ylabel('Rate (Hz)')
    ss4.Box = 'off';
%     ss4.YAxis.Color = ColKDE(1,:);
    
    
    
%     figure()
%     [z,p,k] = butter(6,60/(FS/2),'low');
%     sos_lowY = zp2sos(z,p,k);
%     [z,p,k] = butter(6,60/(1000/2),'low');
%     sos_lowX = zp2sos(z,p,k);
%     YPerStimLow = YPerStim;
%     XAmp_MatLow = XAmp_Mat;
%     for Stim=1:length(YPerStim)
%         YPerStimLow{Stim} = filtfilt(sos_lowY,1,YPerStim{Stim});
%         XAmp_MatLow(Stim,:) = filtfilt(sos_lowX,1,XAmp_Mat(Stim,:));
%     end
%     yyaxis left
%     cla
%     AvAmp = [zeros(1, Delay(1)) mean(XAmp_MatLow,1) zeros(1, Delay(2))];
%     Ampx = -(Delay(1)):1:max(Duration(Indices))+Delay(2)-1;
%     StdAmp = [zeros(1, Delay(1)) std(XAmp_MatLow,0,1)./(size(XAmp_MatLow,1))^0.5 zeros(1, Delay(2))];
%     shadedErrorBar(Ampx,AvAmp,StdAmp,{'-','Color',[0.6350, 0.0780, 0.1840, 0.7], 'LineWidth',2})
%     ylabel('Vocalization Amplitude')
%     
%     yyaxis right
%     cla
%     shadedErrorBar(YPerStimt{1}, mean(cell2mat(YPerStimLow'),1).*1000, 1000*std(cell2mat(YPerStimLow'),0,1)./(length(YPerStimLow)).^0.5, {'-','Color',ColKDE, 'LineWidth',2})
% %     plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
% %     hold on
% %     shadedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
%     hold on
%     VL = vline(0,':k');
%     VL.LineWidth = 2;
%     hold off
%     xlim(XLIM)
%     xlabel('Time (ms)')
%     ylabel('Rate (Hz)')
%     
%     % Calculate the Amplitude MRF
%     AmpWindow = [0 50];
%     NStims=length(Indices);
%     AllAmp = [zeros(NStims, Delay(1)) XAmp_MatLow zeros(NStims, Delay(2))];
%     YPerStimLowMat = cell2mat(YPerStimLow');
%     StepIndices = find((YPerStimt{1}>(-Delay(1)-AmpWindow(1))).*(YPerStimt{1}<(max(YPerStimt{1})-AmpWindow(2))));
%     Nsteps = length(StepIndices);
%     ConvMat = nan(Nsteps*NStims,length(AmpWindow(1):AmpWindow(2)));
%     
%     
%     for oo=1:NStims
%         for step = 1:Nsteps
%             XAmp_startInd = YPerStimt{1}(StepIndices(step)) + Delay(1);
%             XAmp_local = AllAmp(oo, XAmp_startInd+(AmpWindow(1):AmpWindow(2)));
%             ConvMat((oo-1)*Nsteps+step,:)= XAmp_local.* YPerStimLow{oo}(StepIndices(step));
%         end
%     end
%     figure()
%     shadedErrorBar(AmpWindow(1) : AmpWindow(2),mean(ConvMat),std(ConvMat,0,1)./(size(ConvMat,1))^0.5,{'-','Color',[0.6350, 0.0780, 0.1840, 0.7], 'LineWidth',2})
end


    function [YPerStim, YPerStimt,FS] = get_y_4Coherence(SAT, Duration,Delay,TR,Overlap)
        DebugFig = 0;
        if nargin<5
            Overlap = 0;
        end
        if length(Delay)==1
            Delay = [Delay Delay];
        end
        % Calculate the time varying rate applying a gaussian window TR on the
        % spike pattern. The spike pattern considered starts -Delay ms
        % before the onset of the vocalization and stops Delay ms after the
        % offset of the vocalization
        YPerStim = cell(1,length(Duration));
        YPerStimt = cell(1,length(Duration));
        % Gaussian window of 2*std equal to TR (68% of Gaussian centered in TR)
        nStd =(max(Duration) + Delay(1) + Delay(2))/10; % before set as 4
        Tau = (TR/2);
        T_pts = (0:2*nStd*Tau) - nStd*Tau; % centered tpoints around the mean = 0 and take data into account up to nstd away on each side
        Expwav = exp(-0.5*(T_pts).^2./Tau^2)/(Tau*(2*pi)^0.5);
        Expwav = Expwav./sum(Expwav);
        % Frequency at which the neural data should be sampled
        FS = round(1/((TR-Overlap).*10^-3));
        % Loop through the stimuli and fill in the matrix
        for stim=1:length(Duration)
            % Time slots for the neural response
            TimeBinsY = -Delay(1) : (Delay(2) + max(Duration));
            SpikePattern = zeros(1,length(TimeBinsY)-1);
            for isp = 1:length(SAT{stim})
                SpikeInd = round(SAT{stim}(isp));
                if (SpikeInd>=-Delay(1)) && (SpikeInd<(Delay(2) + max(Duration)))
                    SpikePattern(SpikeInd + Delay(1) +1) = SpikePattern(SpikeInd + Delay(1) +1) +1;
                end
            end
            
            % Convolve with Gaussian to obtain our smooth time varying spike train
            % and resample if necessary
            if FS == 1000
                YPerStim{stim} = conv(SpikePattern, Expwav,'same');
                YPerStimt{stim} = TimeBinsY;
            else
                YPerStim_local = conv(SpikePattern, Expwav,'same');
                if sum(YPerStim_local)>0
                    YPerStim_local = YPerStim_local/sum(YPerStim_local)*sum(SpikePattern); % Make sure we keep the right number of sipkes after convolution!
                end
                
                % resampling function is really doing weird things at edges...
                % doing my own resampling
                TimeBinsYOnsetInd = 1 :(TR-Overlap): (Delay(2) + Delay(1) + max(Duration)); % These are slightly different than in get_Y_4GLM, because the times slot are used as indices in the vector and not as actuel time values!
                TimeBinsYOffsetInd = TimeBinsYOnsetInd + TR -1;
                TimeBinsYOnsetInd = TimeBinsYOnsetInd(TimeBinsYOffsetInd<=(Delay(2) + Delay(1) + max(Duration))); % Only keep windows that are within the call
                TimeBinsYOffsetInd = TimeBinsYOffsetInd(TimeBinsYOffsetInd<=(Delay(2) + Delay(1) + max(Duration))); % Only keep windows that are within the call
                
                
                
                TimeBinsYOnset = -Delay(1) :(TR-Overlap): (Delay(2) + max(Duration));
                TimeBinsYOffset = TimeBinsYOnset + TR;
                TimeBinsYOnset = TimeBinsYOnset(TimeBinsYOffset<=(Delay(2) + max(Duration))); % Only keep windows that are within the call
                TimeBinsYOffset = TimeBinsYOffset(TimeBinsYOffset<=(Delay(2) + max(Duration))); % Only keep windows that are within the call
                YPerStimt{stim} = TimeBinsYOnset + (TimeBinsYOffset - TimeBinsYOnset)/2;
                
                YPerStim_resamp = nan(TR, length(TimeBinsYOnsetInd));
                for tt=1:length(TimeBinsYOnsetInd)
                    YPerStim_resamp(:,tt) = YPerStim_local(TimeBinsYOnsetInd(tt): TimeBinsYOffsetInd(tt))';
                end
                YPerStim{stim} = mean(YPerStim_resamp);
                
                
                if DebugFig
                    figure(200) %#ok<UNRCH>
                    clf
                    plot((-Delay(1)+0.5):(Duration(stim)+Delay(2)),YPerStim_local, 'LineWidth',2)
                    xlabel('Time ms')
                    ylabel('Spike Rate mHz (/ms)')
                    title(sprintf('Stim %d/%d',stim,length(Duration)));
                    %         RemainTime = length(XPerStim_temp) - (length(XPerStim{stim})-1)*(1/Fs*10^3);
                    %         XPerStimt{stim} = RemainTime/2+(1/Fs*10^3)*(0:(length(XPerStim{stim})-1));
                    hold on
                    plot(YPerStimt{stim},YPerStim{stim}, 'LineWidth',2)
                    legend({'original' 'resampled'}, 'AutoUpdate','off')
                    hold on
                    SpikeTimes = TimeBinsY(logical(SpikePattern));
                    for ss = 1:length(SpikeTimes)
                        V=vline(SpikeTimes(ss), 'k-');
                        V.LineWidth = 2;
                        hold on
                    end
                    pause(1)
                end
                
            end
            
            %     % change zero values for the smallest value under matlab.
            %     if sum(YPerStim{stim}==0)
            %         MinData = min(YPerStim{stim}(YPerStim{stim} ~=0));
            %         if ~isempty(MinData)
            %             YPerStim{stim}(YPerStim{stim}==0)=min(MinData,realmin('double'));
            %         else
            %             YPerStim{stim}(YPerStim{stim}==0)=realmin('double');
            %         end
            %     end
            
            % Make sure that the output mean(Y) = input mean(Y)
            %     if abs(round(sum(YPerStim{stim})*TR) - sum(SpikePattern))>TR/5
            %         warning('discrepancy in spike rate calculations larger than TR/5= %d?', TR/5)
            %         keyboard
            %     end
            
            %     if any(YPerStim{stim}<0)
            %         keyboard
            %     end
            
            
            %         % Time slots for the neural response
            %         TimeBinsY = -(Delay) : TR: (Delay + Duration(stim));
            %         YPerStim{stim} = nan(1,length(TimeBinsY)-1);
            %         for tt=1:(length(TimeBinsY)-1)
            %             % Find the number of spikes
            %             YPerStim{stim}(tt) = sum( (SAT{stim}>=TimeBinsY(tt)) .* (SAT{stim}<TimeBinsY(tt+1)));
            %         end
        end
        
    end

end
