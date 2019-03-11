function [] = plot_psth_voc_and_behav(SpikeTrainsBehav, SpikeTrainsVoc, Loggers_dir,Date, NeuroLoggerID, Flags, Delay, KDE_Cal)
% AudioLoggerID = ID of the audio logger that the targeting animal is
% wearing
% NeuroLoggerID = ID of the neural logger that the targeting animal is
% wearing
% Flags = whether to print PSTH of Tetrode (Flags(1)=1) and/or Single units
% (Flags(2)=1))
if nargin<7
    KDE_Cal = 0;
end

% Reorganizing input Data
if Flags(1)
    SpikesTTimes_VocCall = SpikeTrainsVoc.SpikesTTimes_VocCall;
    SpikesTTimes_HearCall = SpikeTrainsVoc.SpikesTTimes_HearCall;
    SpikesTTimes_Behav = SpikeTrainsBehav.SpikesTTimes_Behav;
    NT = size(SpikesTTimes_Behav{1},2);
    if KDE_Cal
        Sum_Psth_KDEfiltered_TVocCall = SpikeTrainsVoc.Sum_Psth_KDEfiltered_TVocCall;
        Sum_Psth_KDEfiltered_THearCall = SpikeTrainsVoc.Sum_Psth_KDEfiltered_THearCall;
        Sum_Psth_KDEfiltered_TBehav = SpikeTrainsBehav.Sum_Psth_KDEfiltered_TBehav;
    end
end
if Flags(2)
    SpikesTimes_VocCall = SpikeTrainsVoc.SpikesTimes_VocCall;
    SpikesTimes_HearCall = SpikeTrainsVoc.SpikesTimes_HearCall;
    SpikesTimes_Behav = SpikeTrainsBehav.SpikesTimes_Behav;
    NU = size(SpikesTimes_Behav{1},2);
    if KDE_Cal
        Sum_Psth_KDEfiltered_VocCall = SpikeTrainsVoc.Sum_Psth_KDEfiltered_VocCall;
        Sum_Psth_KDEfiltered_HearCall = SpikeTrainsVoc.Sum_Psth_KDEfiltered_HearCall;
        Sum_Psth_KDEfiltered_Behav = SpikeTrainsBehav.Sum_Psth_KDEfiltered_Behav;
    end
end
VocDuration = SpikeTrainsVoc.VocDuration;
HearDuration = SpikeTrainsVoc.HearDuration;
UActionBehav = SpikeTrainsBehav.UActionBehav;
HearOnlyInd = SpikeTrainsVoc.HearOnlyInd;

MaxInstances = 150; % Max number of instances to plot for behavioral events other than hearing and vocalizing
PSTH_spacer=10;
BehaviorLegendX = -Delay - 50;
%% Now plotting PSTH
fprintf(1, 'Plotting PSTH\n')
% Now plot Raster for tetrodes
ColorCode = get(groot,'DefaultAxesColorOrder');
if Flags(1)
    for uu=1:NT
        fprintf(1, 'Tetrode %d/%d...\n',uu,NT)
        Fig = figure();
        if KDE_Cal
            s2=subplot(3,1,2);
            p2 = get(s2,'Position');
            s2.XTick = [];
            s2.XTickLabel = [];
            s2.YTick = [];
            s2.YTickLabel = [];
            s1=subplot(3,1,1);
            p1 = get(s1,'Position');
            p1(1) = p2(1);
            p1(2) = p2(2);
            p1(4) = p1(4) + p2(4); 
            set(s1,'pos',p1);
        end
        RowCount = 0;
        
        % plotting behavioral action in freely interacting bats
        for bb=1:length(UActionBehav)
            fprintf(1,' -> Raster %s\n', UActionBehav{bb})
            NBehav = min(size(SpikesTTimes_Behav{bb},1), MaxInstances);
            
            % Plotting spikes
            for cc=1:NBehav
                hold on
                % first plot a colored line to code the behavior
                plot([-Delay 3*Delay], RowCount+cc-[0.5 0.5], '-','LineWidth',250/(NBehav*length(UActionBehav)),'Color', [ColorCode(2+bb,:) 0.1])
                % then the spikes
                for spike=1:length(SpikesTTimes_Behav{bb}{cc,uu})
                    hold on
                    plot(SpikesTTimes_Behav{bb}{cc,uu}(spike)*ones(2,1)-Delay, RowCount+cc-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
                %                 yyaxis right
                %                 plot(Psth_KDEfiltered_TVocCall_t{cc,uu}, Psth_KDEfiltered_TVocCall{cc,uu}/max(Psth_KDEfiltered_TVocCall_scalef(:,uu))+cc-1, 'r-', 'LineWidth',2)
            end
            % plotting the legend for that particular behavior
            h=text(BehaviorLegendX, RowCount + NBehav/3 ,sprintf('%s', UActionBehav{bb}));
            set(h,'Rotation',90);
            RowCount = RowCount + NBehav + PSTH_spacer; % Increment the row count and add PSTH_spacer to leave PSTH_spacer row space between 2 behavior types
            
        end
        
        % plotting hearing data during conditioning
        if isnan(HearOnlyInd)
            HearOnlyDuration = Inf;
            fprintf(1,' ->No hearing data from the operant conditioning\n')
        else
            fprintf(1,' -> Raster hearing\n')
            HearOnlyDuration=HearDuration(HearOnlyInd); % Only taking heard vocalizations when the bat is not vocalizing itself.
            % We want to plot PSTH with increasing duration of vocalizations
            [~, IDurH] = sort(HearOnlyDuration, 'descend');
            for hh=1:length(HearOnlyInd)
                cc = HearOnlyInd(IDurH(hh));
%                 yyaxis left
                hold on
                plot([0 HearDuration(cc)], RowCount+hh-[0.5 0.5], '-','LineWidth',250/(length(HearOnlyInd)+ RowCount),'Color', [0.8 0.8 1])
                %         for dd=1:size(DataDeletion_HearCall{cc},1)
                %             hold on
                %             plot(DataDeletion_HearCall{cc}(dd,:), hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 0.8]) % RF artefact period
                %         end
                for spike=1:length(SpikesTTimes_HearCall{cc,uu})
                    hold on
                    plot(SpikesTTimes_HearCall{cc,uu}(spike)*ones(2,1), RowCount+hh-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
                %             yyaxis right
                %             plot(Psth_KDEfiltered_THearCall_t{cc,uu}, Psth_KDEfiltered_THearCall{cc,uu}/max(Psth_KDEfiltered_THearCall_scalef(HearOnlyInd,uu))+hh-1, 'r-', 'LineWidth',2)
            end
            % plotting a line at onset
            line([0 0], [RowCount-0.5 RowCount+length(HearOnlyInd)+0.5], 'color','b','LineStyle','-', 'LineWidth',1)
            % plotting the legend for that particular behavior
            h=text(BehaviorLegendX, RowCount + length(HearOnlyInd)/3 ,'Hearing');
            set(h,'Rotation',90);
            RowCount = RowCount + length(HearOnlyInd) + PSTH_spacer; % Increment the row count and add 1 to leave one row space between 2 behavior types
        end
        
        % plotting vocalization production data during conditioning
        % We want to plot PSTH with increasing duration of vocalizations
        fprintf(1,' -> Raster vocalizing\n')
        [~, IDurV] = sort(VocDuration, 'descend');
        for hh=1:length(VocDuration)
            cc = IDurV(hh);
%             yyaxis left
            hold on
            plot([0 VocDuration(cc)], RowCount+hh-[0.5 0.5], '-','LineWidth',250/(length(VocDuration)+ RowCount),'Color', [1 0.8 0.8])
            %         for dd=1:size(DataDeletion_HearCall{cc},1)
            %             hold on
            %             plot(DataDeletion_HearCall{cc}(dd,:), hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 0.8]) % RF artefact period
            %         end
            for spike=1:length(SpikesTTimes_VocCall{cc,uu})
                hold on
                plot(SpikesTTimes_VocCall{cc,uu}(spike)*ones(2,1), RowCount+hh-[0.9 0.1], 'k-', 'LineWidth',1)
            end
            hold on
            %             yyaxis right
            %             plot(Psth_KDEfiltered_THearCall_t{cc,uu}, Psth_KDEfiltered_THearCall{cc,uu}/max(Psth_KDEfiltered_THearCall_scalef(HearOnlyInd,uu))+hh-1, 'r-', 'LineWidth',2)
        end
        % plotting a line at onset
        line([0 0], [RowCount-0.5 RowCount+length(VocDuration)+0.5], 'Color','r','LineStyle','-', 'LineWidth',1.5)
        % plotting the legend for that particular behavior
        h=text(BehaviorLegendX, RowCount + length(VocDuration)/3 ,'Vocalizing');
        set(h,'Rotation',90);
        RowCount = RowCount + length(VocDuration)+1;
        
        % Set the parameters of the figure
        XLIM = [-Delay+10 min(min(HearOnlyDuration), min(VocDuration))+Delay-10];
        xlabel('Time (ms)')
        ylim([0 RowCount])
        xlim(XLIM)
        title(sprintf('Raster %s Tetrode %d on %s', NeuroLoggerID, uu, Date))
        hold off
        
        % Plot KDE if requested
        if KDE_Cal
            MaxSR = nan(length(UActionBehav)+2,1);
            subplot(3,1,3)
            % First plot lines to get the legend right
            for bb=1:length(UActionBehav)
                hold on
                plot(Sum_Psth_KDEfiltered_TBehav{uu,bb}{1}-Delay, Sum_Psth_KDEfiltered_TBehav{uu,bb}{2},'Color',ColorCode(2+bb,:), 'LineStyle','-', 'LineWidth',2)
                Ind = ((Sum_Psth_KDEfiltered_TBehav{uu,bb}{1}-Delay)>XLIM(1)) .* ((Sum_Psth_KDEfiltered_TBehav{uu,bb}{1}-Delay) < XLIM(2));
                if ~isempty(Ind)
                    MaxSR(bb) = max(Sum_Psth_KDEfiltered_TBehav{uu,bb}{2}(logical(Ind)));
                end
            end
            
            if (length(HearOnlyInd)>1) && KDE_Cal
                hold on
                plot(Sum_Psth_KDEfiltered_THearCall{uu,1}, Sum_Psth_KDEfiltered_THearCall{uu,2}, 'b-', 'LineWidth',2)
                Ind = ((Sum_Psth_KDEfiltered_THearCall{uu,1}-Delay)>XLIM(1)) .* ((Sum_Psth_KDEfiltered_THearCall{uu,1}-Delay) < XLIM(2));
                MaxSR(bb+1) = max(Sum_Psth_KDEfiltered_THearCall{uu,2}(logical(Ind)));
            end
            hold on
            plot(Sum_Psth_KDEfiltered_TVocCall{uu,1}, Sum_Psth_KDEfiltered_TVocCall{uu,2}, 'r-', 'LineWidth',2)
            Ind = ((Sum_Psth_KDEfiltered_TVocCall{uu,1}-Delay)>XLIM(1)) .* ((Sum_Psth_KDEfiltered_TVocCall{uu,1}-Delay) < XLIM(2));
            MaxSR(bb+2) = max(Sum_Psth_KDEfiltered_TVocCall{uu,2}(logical(Ind)));
            xlabel('Time (ms)')
            ylabel('Spike rate (/ms)')
            % Then plots the lines with error bars
            if (length(HearOnlyInd)>1) && KDE_Cal
                legend([UActionBehav; 'Hearing'; 'Vocalizing'], 'AutoUpdate', 'off','location','southoutside','Orientation','horizontal')
            else
                legend([UActionBehav; 'Vocalizing'], 'AutoUpdate', 'off','location','southoutside','Orientation','horizontal')
            end
            for bb=1:length(UActionBehav)
                fprintf(1,' -> KDE %s\n', UActionBehav{bb})
                hold on
                shadedErrorBar(Sum_Psth_KDEfiltered_TBehav{uu,bb}{1}-Delay, Sum_Psth_KDEfiltered_TBehav{uu,bb}{2}, Sum_Psth_KDEfiltered_TBehav{uu,bb}{3}, {'Color',ColorCode(2+bb,:), 'LineStyle','-', 'LineWidth',2})
            end
            
            if (length(HearOnlyInd)>1) && KDE_Cal
                hold on
                fprintf(1, ' -> KDE Hearing\n')
                shadedErrorBar(Sum_Psth_KDEfiltered_THearCall{uu,1}, Sum_Psth_KDEfiltered_THearCall{uu,2}, Sum_Psth_KDEfiltered_THearCall{uu,3}, {'b-', 'LineWidth',2})
            end
            hold on
            fprintf(1, ' -> KDE Vocalizing\n')
            shadedErrorBar(Sum_Psth_KDEfiltered_TVocCall{uu,1}, Sum_Psth_KDEfiltered_TVocCall{uu,2}, Sum_Psth_KDEfiltered_TVocCall{uu,3}, {'r-', 'LineWidth',2})
            YLIM = get(gca,'YLim');
            YLIM(1) = 0;
            YLIM(2) = max(MaxSR)*1.1;
            ylim(YLIM);
            xlim(XLIM)
            hold off
        end
        
        % Save figure
        orient(Fig,'portrait')
        Fig.PaperPositionMode = 'auto';
        set(Fig,'PaperOrientation','portrait');
        %             set(Fig,'PaperUnits','normalized');
        %             set(Fig,'PaperPosition', [50 50 1200 800]);
        %             pause()
        if KDE_Cal
            print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_ALL_PSTH_KDE_Tetrode%d.pdf', Date, NeuroLoggerID,uu)),'-dpdf','-fillpage')
        else
            print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_ALL_PSTH_Tetrode%d.pdf', Date, NeuroLoggerID,uu)),'-dpdf','-fillpage')
        end
        close all
        
    end
    
end

        
 %% Now plot Raster for single units
if Flags(2)
    for uu=1:NU
        fprintf(1, 'Single Unit %d/%d...\n',uu,NU)
        Fig = figure();
        if KDE_Cal
            s2=subplot(3,1,2);
            p2 = get(s2,'Position');
            s2.XTick = [];
            s2.XTickLabel = [];
            s2.YTick = [];
            s2.YTickLabel = [];
            s1=subplot(3,1,1);
            p1 = get(s1,'Position');
            p1(1) = p2(1);
            p1(2) = p2(2);
            p1(4) = p1(4) + p2(4); 
            set(s1,'pos',p1);
        end
        RowCount = 0;
        
        % plotting behavioral action in freely interacting bats
        for bb=1:length(UActionBehav)
            NBehav = min(size(SpikesTimes_Behav{bb},1), MaxInstances);
            fprintf(1,' -> %s\n', UActionBehav{bb})
            % Plotting spikes
            for cc=1:NBehav
                hold on
                % First plot a colored line to indicate the behavior
                plot([-Delay 3*Delay], RowCount+cc-[0.5 0.5], '-','LineWidth',250/(NBehav*length(UActionBehav)),'Color', [ColorCode(2+bb,:) 0.1])
                hold on
                % then the spikes
                for spike=1:length(SpikesTimes_Behav{bb}{cc,uu})
                    hold on
                    plot(SpikesTimes_Behav{bb}{cc,uu}(spike)*ones(2,1)-Delay, RowCount+cc-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
                %                 yyaxis right
                %                 plot(Psth_KDEfiltered_TVocCall_t{cc,uu}, Psth_KDEfiltered_TVocCall{cc,uu}/max(Psth_KDEfiltered_TVocCall_scalef(:,uu))+cc-1, 'r-', 'LineWidth',2)
            end
            % plotting the legend for that particular behavior
            h=text(BehaviorLegendX, RowCount + NBehav/3 ,sprintf('%s', UActionBehav{bb}));
            set(h,'Rotation',90);
            RowCount = RowCount + NBehav + PSTH_spacer; % Increment the row count and add PSTH_spacer to leave PSTH_spacer row space between 2 behavior types
            
        end
        
        % plotting hearing data during conditioning
        if isnan(HearOnlyInd)
            fprintf(1,' -> No hearing data from the operant conditioning\n')
        else
            fprintf(1,' -> hearing\n')
            HearOnlyDuration=HearDuration(HearOnlyInd); % Only taking heard vocalizations when the bat is not vocalizing itself.
            % We want to plot PSTH with increasing duration of vocalizations
            [~, IDurH] = sort(HearOnlyDuration, 'descend');
            for hh=1:length(HearOnlyInd)
                cc = HearOnlyInd(IDurH(hh));
%                 yyaxis left
                hold on
                plot([0 HearDuration(cc)], RowCount+hh-[0.5 0.5], '-','LineWidth',250/(length(HearOnlyInd)+ RowCount),'Color', [0.8 0.8 1])
                %         for dd=1:size(DataDeletion_HearCall{cc},1)
                %             hold on
                %             plot(DataDeletion_HearCall{cc}(dd,:), hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 0.8]) % RF artefact period
                %         end
                for spike=1:length(SpikesTimes_HearCall{cc,uu})
                    hold on
                    plot(SpikesTimes_HearCall{cc,uu}(spike)*ones(2,1), RowCount+hh-[0.9 0.1], 'k-', 'LineWidth',1)
                end
                hold on
                %             yyaxis right
                %             plot(Psth_KDEfiltered_THearCall_t{cc,uu}, Psth_KDEfiltered_THearCall{cc,uu}/max(Psth_KDEfiltered_THearCall_scalef(HearOnlyInd,uu))+hh-1, 'r-', 'LineWidth',2)
            end
            % plotting a line at onset
            line([0 0], [RowCount-0.5 RowCount+length(HearOnlyInd)+0.5],'color','b','LineStyle','-', 'LineWidth',1)
            % plotting the legend for that particular behavior
            h=text(BehaviorLegendX, RowCount + length(HearOnlyInd)/3 ,'Hearing');
            set(h,'Rotation',90);
            RowCount = RowCount + length(HearOnlyInd) + PSTH_spacer; % Increment the row count and add PSTH_spacer to leave PSTH_spacer row space between 2 behavior types
        end
        
        % plotting vocalization production data during conditioning
        % We want to plot PSTH with increasing duration of vocalizations
        fprintf(1,' -> vocalizing\n')
        [~, IDurV] = sort(VocDuration, 'descend');
        for hh=1:length(VocDuration)
            cc = IDurV(hh);
%             yyaxis left
            hold on
            plot([0 VocDuration(cc)], RowCount+hh-[0.5 0.5], '-','LineWidth',250/(length(VocDuration)+ RowCount),'Color', [1 0.8 0.8])
            %         for dd=1:size(DataDeletion_HearCall{cc},1)
            %             hold on
            %             plot(DataDeletion_HearCall{cc}(dd,:), hh-[0.5 0.5], '-','LineWidth',250/length(HearOnlyInd),'Color', [0.8 0.8 0.8]) % RF artefact period
            %         end
            for spike=1:length(SpikesTimes_VocCall{cc,uu})
                hold on
                plot(SpikesTimes_VocCall{cc,uu}(spike)*ones(2,1), RowCount+hh-[0.9 0.1], 'k-', 'LineWidth',1)
            end
            hold on
            %             yyaxis right
            %             plot(Psth_KDEfiltered_THearCall_t{cc,uu}, Psth_KDEfiltered_THearCall{cc,uu}/max(Psth_KDEfiltered_THearCall_scalef(HearOnlyInd,uu))+hh-1, 'r-', 'LineWidth',2)
        end
        % plotting a line at onset
        line([0 0], [RowCount-0.5 RowCount+length(VocDuration)+0.5], 'color','r','LineStyle','-', 'LineWidth',1)
        % plotting the legend for that particular behavior
        h=text(BehaviorLegendX, RowCount + length(VocDuration)/3 ,'Vocalizing');
        set(h,'Rotation',90);
        RowCount = RowCount + length(VocDuration) + 1;
        
        % Set the parameters of the figure
        XLIM = [-Delay+10 min(min(HearOnlyDuration), min(VocDuration))+Delay-10];
        xlabel('Time (ms)')
        ylim([0 RowCount])
        xlim(XLIM)
        title(sprintf('Raster %s Single Unit %d on %s', NeuroLoggerID, uu, Date))
        hold off
        
        % Plot KDE if requested 
        if KDE_Cal
            MaxSR = nan(length(UActionBehav)+2,1);
            subplot(3,1,3)
            % First plot the lines to get the legend right
            for bb=1:length(UActionBehav)
                hold on
                plot(Sum_Psth_KDEfiltered_Behav{uu,bb}{1}-Delay, Sum_Psth_KDEfiltered_Behav{uu,bb}{2}, 'Color',ColorCode(2+bb,:), 'LineStyle','-', 'LineWidth',2)
                Ind = ((Sum_Psth_KDEfiltered_Behav{uu,bb}{1}-Delay)>XLIM(1)) .* ((Sum_Psth_KDEfiltered_Behav{uu,bb}{1}-Delay) < XLIM(2));
                if ~isempty(Ind)
                    MaxSR(bb) = max(Sum_Psth_KDEfiltered_Behav{uu,bb}{2}(logical(Ind)));
                end
            end
            
            if (length(HearOnlyInd)>1) && KDE_Cal
                hold on
                plot(Sum_Psth_KDEfiltered_HearCall{uu,1}, Sum_Psth_KDEfiltered_HearCall{uu,2}, 'b-', 'LineWidth',2)
                Ind = ((Sum_Psth_KDEfiltered_HearCall{uu,1}-Delay)>XLIM(1)) .* ((Sum_Psth_KDEfiltered_HearCall{uu,1}-Delay) < XLIM(2));
                MaxSR(bb+1) = max(Sum_Psth_KDEfiltered_HearCall{uu,2}(logical(Ind)));
            end
            hold on
            plot(Sum_Psth_KDEfiltered_VocCall{uu,1}, Sum_Psth_KDEfiltered_VocCall{uu,2}, 'r-', 'LineWidth',2)
            Ind = ((Sum_Psth_KDEfiltered_VocCall{uu,1}-Delay)>XLIM(1)) .* ((Sum_Psth_KDEfiltered_VocCall{uu,1}-Delay) < XLIM(2));
            MaxSR(bb+2) = max(Sum_Psth_KDEfiltered_VocCall{uu,2}(logical(Ind)));
            xlabel('Time (ms)')
            ylabel('Spike rate (/ms)')
            if (length(HearOnlyInd)>1) && KDE_Cal
                legend([UActionBehav; 'Hearing'; 'Vocalizing'], 'AutoUpdate','off','location','southoutside','Orientation','horizontal')
            else
                legend([UActionBehav; 'Vocalizing'], 'AutoUpdate','off','location','southoutside','Orientation','horizontal')
            end
            for bb=1:length(UActionBehav)
                fprintf(1,' -> KDE %s\n', UActionBehav{bb})
                hold on
                shadedErrorBar(Sum_Psth_KDEfiltered_Behav{uu,bb}{1}-Delay, Sum_Psth_KDEfiltered_Behav{uu,bb}{2}, Sum_Psth_KDEfiltered_Behav{uu,bb}{3}, {'Color',ColorCode(2+bb,:), 'LineStyle','-', 'LineWidth',2})
            end
            
            if (length(HearOnlyInd)>1) && KDE_Cal
                hold on
                fprintf(1,' -> KDE Hearing\n')
                shadedErrorBar(Sum_Psth_KDEfiltered_HearCall{uu,1}, Sum_Psth_KDEfiltered_HearCall{uu,2}, Sum_Psth_KDEfiltered_HearCall{uu,3}, {'b-', 'LineWidth',2})
            end
            hold on
            fprintf(1,' -> KDE vocalizing\n')
            shadedErrorBar(Sum_Psth_KDEfiltered_VocCall{uu,1}, Sum_Psth_KDEfiltered_VocCall{uu,2}, Sum_Psth_KDEfiltered_VocCall{uu,3}, {'r-', 'LineWidth',2})
            YLIM = get(gca,'YLim');
            YLIM(1) = 0;
            YLIM(2) = max(MaxSR)*1.1;
            ylim(YLIM);
            xlim(XLIM)
            hold off
        end
        
        % Save figure
        orient(Fig,'portrait')
        Fig.PaperPositionMode = 'auto';
        set(Fig,'PaperOrientation','portrait');
        %             set(Fig,'PaperUnits','normalized');
        %             set(Fig,'PaperPosition', [50 50 1200 800]);
        %             pause()
        if KDE_Cal
            print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_ALL_PSTH_KDE_SU%d.pdf', Date, NeuroLoggerID,uu)),'-dpdf','-fillpage')
        else
            print(Fig,fullfile(Loggers_dir,sprintf('%s_%s_ALL_PSTH_SU%d.pdf', Date, NeuroLoggerID,uu)),'-dpdf','-fillpage')
        end
        close all
        
    end
    
end
fprintf(1, 'DONE\n')

end
