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
IndVocF = find(contains(Data.ExpType, 'F'));
IndVocPD = intersect(IndVocP, IndVocD);
IndVocPDO = intersect(IndVocPD, IndVocO);
IndVocPDF = intersect(IndVocPD, IndVocF);
IndVocHD = intersect(IndVocH, IndVocD);
IndVocHDO = intersect(IndVocHD, IndVocO);
IndVocHDF = intersect(IndVocHD, IndVocF);
%% Time Raster plot alligned to vocalization production onset/offset self vocalizations Operant + Free First voc of sequence only
if ~isempty(IndVocPD) && ~isempty(IndVocPDO) && ~isempty(IndVocPDF)
    Fig1 = figure();% TrCol = [0.9290, 0.6940, 0.1250];BaCol = [1, 0, 0];
    Color = [0/255 191/255 255/255].*contains(Data.What, 'Tr') + [1 0.7 0.7].*contains(Data.What, 'Ba');
    if isfield(Data.KDE_onset,'SelfVocAll')
        ColKDE = [186/255 85/255 211/255];
        timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPD,Color,Data.KDE_onset.SelfVocAll,Data.KDE_offset.SelfVocAll,ColKDE,Data.RewardTime,DurOrd);
    else
        timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPD,Color);
    end
    suplabel(sprintf('CALLS FROM SUBJECT O and F   %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
    print(Fig1,fullfile(OutputPath,sprintf('%s_RasterVocSelf_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage');
end

%% Time Raster plot alligned to vocalization perception onset/offset others vocalizations Operant + Free First voc of sequence only
if ~isempty(IndVocHD) && ~isempty(IndVocHDO) && ~isempty(IndVocHDF)
    Fig2 = figure();
    Color = [0 0.3 0].*contains(Data.What, 'Tr') + repmat([0.7 0.7 1], length(Data.What),1);
    if isfield(Data.KDE_onset, 'OthersVocAll')
        timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHD,Color, Data.KDE_onset.OthersVocAll,Data.KDE_offset.OthersVocAll,[0.4 0.45 1],Data.RewardTime,DurOrd)
    else
        timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHD,Color)
    end
    suplabel(sprintf('CALLS FROM OTHERS O and F    %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
    print(Fig2,fullfile(OutputPath,sprintf('%s_RasterVocOthers_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
end

%% Time Raster plot alligned to vocalization production onset/offset during Operant conditioning First voc of sequence only
if ~isempty(IndVocPDO)
    Fig6 = figure();
    Color = [0/255 191/255 255/255].*contains(Data.What, 'Tr') + [1 0.7 0.7].*contains(Data.What, 'Ba');
    if isfield(Data.KDE_onset, 'SelfVocOp')
        ColKDE = [186/255 85/255 211/255];
        timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPDO,Color,Data.KDE_onset.SelfVocOp,Data.KDE_offset.SelfVocOp,ColKDE,Data.RewardTime,DurOrd)
    else
        timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPDO,Color)
    end
    suplabel(sprintf('CALLS FROM SUBJECT OPERANT   %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
    print(Fig6,fullfile(OutputPath,sprintf('%s_RasterVocSelfOp_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
end

%% Time Raster plot alligned to vocalization perception onset/offset during operant conditioning First voc of sequence only
if ~isempty(IndVocHDO)
    Fig7 = figure();
    Color = [0 0.3 0].*contains(Data.What, 'Tr') + repmat([0.7 0.7 1], length(Data.What),1);
    if isfield(Data.KDE_onset,'OthersVocOp')
        timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHDO,Color,Data.KDE_onset.OthersVocOp,Data.KDE_offset.OthersVocOp,[0.4 0.45 1],Data.RewardTime,DurOrd)
    else
        timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHDO,Color)
    end
    suplabel(sprintf('CALLS FROM OTHERS OPERANT    %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
    print(Fig7,fullfile(OutputPath,sprintf('%s_RasterVocOthersOp_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
end

%% Time Raster plot alligned to vocalization production onset/offset during Free session First voc of sequence only
if ~isempty(IndVocPDF)
    Fig8 = figure();
    Color = [0/255 191/255 255/255].*contains(Data.What, 'Tr') + [1 0.7 0.7].*contains(Data.What, 'Ba');
    if isfield(Data.KDE_onset,'SelfVocFr')
        ColKDE = [186/255 85/255 211/255];
        timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPDF,Color,Data.KDE_onset.SelfVocFr,Data.KDE_offset.SelfVocFr,ColKDE)
    else
        timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocPDF,Color)
    end
    suplabel(sprintf('CALLS FROM SUBJECT FREE SESSION    %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
    print(Fig8,fullfile(OutputPath,sprintf('%s_RasterVocSelfFr_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
end

%% Time Raster plot alligned to vocalization perception onset/offset during Free session First voc of sequence only
if ~isempty(IndVocHDF)
    Fig9 = figure();
    Color = [0 0.3 0].*contains(Data.What, 'Tr') + repmat([0.7 0.7 1], length(Data.What),1);
    if isfield(Data.KDE_onset,'OthersVocFr')
        timerasterkde(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHDF,Color, Data.KDE_onset.OthersVocFr,Data.KDE_offset.OthersVocFr,[0.4 0.45 1])
    else
        timeraster(Data.SpikesArrivalTimes_Behav,Data.Duration,Delay,IndVocHDF,Color)
    end
    suplabel(sprintf('CALLS FROM OTHERS FREE SESSION    %s on %s Raster T%s SS%s %s',SubjectID, Date, NeuralInputID{1},NeuralInputID{3},NeuralInputID{2}),'t');
    print(Fig9,fullfile(OutputPath,sprintf('%s_RasterVocOthersFr_%d.pdf', FileNameBase, Delay(1))),'-dpdf','-fillpage')
end

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

%% INTERNAL FUNCTION
    function timeraster(SpikesArrivalTimes,Duration,Delay,Indices, Color)
        % We want to plot data with increasing duration of
        % vocalizations.
        [~,IDur] = sort(Duration(Indices));

        subplot(1,2,1)

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
        hold off
        
        
        % then alligned to vocalization offset
        subplot(1,2,2)

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
        hold off
    end


function timerasterkde(SpikesArrivalTimes,Duration,Delay,Indices, Color, Dat1, Dat2, Col1, RewardTime,DurOrd)
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
    
        subplot(4,2,[1 3 5])

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
        subplot(4,2,[2 4 6])

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
end
