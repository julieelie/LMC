function [AllActions, UActionText]= get_logger_data_event(ExtData_dir, Date, Buffer, LoggerID, Output_dir)
% AllActions is a cell array of length the number of actions containing
% each 2 columns giving the onset and offset times in ms of behavioral
% events
if nargin<3
    Buffer = 500; %ms
end

AllSep = strfind(ExtData_dir, filesep);
if nargin<4
    if length(ExtData_dir)==AllSep(end)
        LoggerID = ExtData_dir((AllSep(end-2)+1):(AllSep(end-1)-1));
    else
        LoggerID = ExtData_dir((AllSep(end-1)+1):(AllSep(end)-1));
    end
end
if nargin<5
    if length(ExtData_dir)==AllSep(end)
        Output_dir = ExtData_dir(1:(AllSep(end-2)));
    else
        Output_dir = ExtData_dir(1:(AllSep(end-1)));
    end
end
fprintf(1, 'Figures and data are saved to %s\n', Output_dir)
MaxEventDur = 1000; %ms
EventfilePath = dir(fullfile(ExtData_dir, sprintf('*%s*EVENTS.mat', Date)));
Eventfile = load(fullfile(EventfilePath.folder, EventfilePath.name));

%% identify all the free text strings reported
FreeTextInd = find(contains(Eventfile.event_types_and_details, 'Free text'));
UFreeText_all = unique(Eventfile.event_types_and_details(FreeTextInd)); %#ok<FNDSB> % unique labels
% Remove 'Free Text' and isolate actions items
SpaceInd = strfind(UFreeText_all, ' ');
StartInd= strfind(UFreeText_all, 'start');
StopInd= strfind(UFreeText_all, 'stop');
UActionText = cell(size(UFreeText_all));
UFullActionText = cell(size(UFreeText_all));
for ii=1:length(SpaceInd)
    UFullActionText{ii} = UFreeText_all{ii}((SpaceInd{ii}(2)+1) :end);
    if ~isempty(StartInd{ii})
        UActionText{ii} = UFreeText_all{ii}((SpaceInd{ii}(2)+1) :(StartInd{ii}-1));
    elseif ~isempty(StopInd{ii})
        UActionText{ii} = UFreeText_all{ii}((SpaceInd{ii}(2)+1) :(StopInd{ii}-1));
    else
        UActionText{ii} = UFreeText_all{ii}((SpaceInd{ii}(2)+1):end);
    end
end
UActionText = unique(UActionText);

%% construct an ethogram of the recording
NAction = length(UActionText);
AllActions = cell(NAction,1);
% Find the reference time to plot the ethogram
RefTime = Eventfile.event_timestamps_usec(find(contains(Eventfile.event_types_and_details, 'Started recording'),1))*10^-3;
% Get the indices of all the starting points
IndStart = find(contains(Eventfile.event_types_and_details, 'start'));
% Get the indices of all the stoping points
IndStop = find(contains(Eventfile.event_types_and_details, 'stop'));
Fig=figure();
ColorCode = get(groot,'DefaultAxesColorOrder');
ColorCode = [ColorCode(1:6,:) ; 0 0 1; ColorCode(7,:); 0.85 0.6940 0.556; 0 0 0; 1 0 0; 0.301 0.078 0.741; ];
% loop through each action
for aa=1:NAction
    fprintf('extracting data of %s\n',UActionText{aa})
    IndAction = find(contains(Eventfile.event_types_and_details, UActionText{aa}));
    IndActionStart = intersect(IndAction, IndStart);
    if isempty(IndActionStart) % This action is ponctual
        fprintf('Ploting ponctual events\n')
        AllActions{aa} = Eventfile.event_timestamps_usec(IndAction)*10^-3;
        hold on
        plot((AllActions{aa}-RefTime)*10^-3, aa, '*','Color', ColorCode(aa,:), 'MarkerSize', 5)
    else
        fprintf('Finding start and end of each event\n')
        IndActionStop = intersect(IndAction, IndStop);
        AllActions{aa} = nan(max(length(IndActionStop), length(IndActionStart)),2);
        % loop through action starts and figure out if a stop is closer to
        % it than the next start
        e_count = 0;
        for ii=1:length(IndActionStart)
            fprintf('potential start %d/%d\n',ii,length(IndActionStart));
            % get the stop closer to the start
            mm = find(min(IndActionStop - IndActionStart(ii)));
            if (IndActionStop(mm) - IndActionStart(ii)) <0
                error('Error in parsing the action, the stop happens before the start\n')
            end
            if ~isempty(IndActionStop) && ((ii==length(IndActionStart)) || (IndActionStop(mm) - IndActionStart(ii)) < diff(IndActionStart(ii:(ii+1))))  % this stop happens before another start, all is good save that pair of event
                e_count= e_count+1;
                AllActions{aa}(e_count,1) = Eventfile.event_timestamps_usec(IndActionStart(ii))*10^-3;
                AllActions{aa}(e_count,2) = Eventfile.event_timestamps_usec(IndActionStop(mm))*10^-3;
                % plot the action
                hold on
                plot((AllActions{aa}(e_count,:)-RefTime)*10^-3, aa*ones(2,1), '-', 'Color', ColorCode(aa,:), 'LineWidth',2);
                if ii<length(IndActionStart)
                    % Find the next start
                    nn = find(min(IndActionStart(ii+1:end) - IndActionStop(mm)));
                    if (IndActionStart(ii+nn) - IndActionStop(mm)) < diff(IndActionStop(mm:(mm+1)))
                        % The next start happens before another stop, all is
                        % good, proceed to the next start and supress that stop
                        IndActionStop(mm) = [];
                    else % There are several stops in a row, get rid of them up to the next start
                        SupInd = find(diff(IndActionStop(mm:end))<(IndActionStart(ii+nn) - IndActionStop(mm)));
                        IndActionStop(mm : (mm+ SupInd)) = [];
                    end
                end
            else 
                % This start misses a stop, do not consider it consider the next start
            end
        end
        AllActions{aa} = AllActions{aa}(1:e_count,:);
    end
end

set(gca,'YLim',[0 NAction+1], 'YTick',1:NAction, 'YTickLabel', UActionText)
xlabel('Time (in s)')
title(sprintf('Ethogram Implanted bat %s', Date))

orient(Fig,'landscape')
Fig.PaperPositionMode = 'auto';
set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
set(Fig,'PaperOrientation','landscape');
set(Fig,'PaperUnits','normalized');
set(Fig,'PaperPosition', [0 0 1 1]);

print(Fig,fullfile(Output_dir,sprintf('%s_%s_Actuogram.pdf', Date(3:end), LoggerID)),'-dpdf')

%% extract the neural activity during the whole recording time
PlayBackOffset = AllActions{contains(UActionText, 'playback')}(end,2);
fprintf('Extract Single Unit spike arrival times for the whole experiment session\n')
Delay = 0;% time in ms to extract data before and after event onset/offset
[~,~, ~, SpikeSU] = extract_timeslot_LFP_spikes(ExtData_dir, [RefTime, PlayBackOffset], Delay, NaN,[0 0 0 1]);
save(fullfile(Output_dir, sprintf('%s_SpikeSU_RecordPlayback.mat', Date(3:end))),'SpikeSU');
Response_samprate = 1;% Sampling rate of the KDE in Hz
Bin_ms = 1000; % size of the KDE binning
t=-Delay: Bin_ms : round((PlayBackOffset-RefTime + Delay)/Bin_ms)*Bin_ms;
Psth_KDEfiltered_SpikeSU = cell(1,length(SpikeSU));
Psth_KDEfiltered_SpikeSU_t = cell(1,length(SpikeSU));
for uu=1:length(SpikeSU)
    fprintf('calculate spike density estimate SU %d/%d\n', uu, length(SpikeSU))
    % calculate the density estimate (Kernel density estimation: L2 risk
    % minimization, with reflection boundary)
    [y,Psth_KDEfiltered_SpikeSU_t{uu},~]=ssvkernel(SpikeSU{uu},t);
%     [Psth_KDEfiltered_SpikeSU_t{uu},y]=Bayes_v1(SpikeSU{uu});
    % y is a density function that sums to 1
    % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 1 s)
    % multiplying by the response sampling rate in Hz gives the expected spike rate in spike/s
    Psth_KDEfiltered_SpikeSU{uu} =  y * length(SpikeSU{uu}) * Response_samprate;
end

Fig=figure();
% Plot KDE
SpikeCol = [1 0.2 0.2 0.3;  0.5 0.7 0.5 0.8; 0.9 0.1 0.2 0.3; [get(groot,'DefaultAxesColorOrder') repmat(0.3,7,1)]];
LegendKDE =cell(length(SpikeSU),1);
for uu=1:length(SpikeSU)
    yyaxis right
    hold on
    plot(Psth_KDEfiltered_SpikeSU_t{uu}*10^-3,zscore(Psth_KDEfiltered_SpikeSU{uu}),'-','Color',SpikeCol(uu,:), 'LineWidth',2, 'DisplayName', sprintf('Single Unit %d',uu))
    LegendKDE{uu}  = sprintf('Single Unit %d',uu);
    hold off    
end
ylabel('Zscored Spike Rate (/s)')
legend(LegendKDE, 'AutoUpdate','off')
YlimRight = get(gca, 'YLim');
set(gca, 'YLim', [0 YlimRight(2)*2])

% Plot actions on top
yyaxis left
for aa=1:NAction
    if size(AllActions{aa},2)==1
        hold on
        plot((AllActions{aa}-RefTime)*10^-3, aa, '*','Color', ColorCode(aa,:), 'MarkerSize', 5)
    elseif size(AllActions{aa},2)==2
        for e_count=1:size(AllActions{aa},1)
            hold on
            plot((AllActions{aa}(e_count,:)-RefTime)*10^-3, aa*ones(2,1), '-', 'Color', ColorCode(aa,:), 'LineWidth',2);
        end
    end
end
xlabel('Time (in s)')
title(sprintf('Ethogram and SU bat %s', Date))

% Plot spikes!
LegendPlot = cell(length(SpikeSU) + length(UActionText),1);
LegendPlot((length(SpikeSU)+1) : end) = UActionText;
yyaxis left
for uu=1:length(SpikeSU)
    hold on
    plot(SpikeSU{uu}*10^-3, ones(length(SpikeSU{uu}),1)*(-uu/(length(SpikeSU)+1)), '.', 'MarkerSize',5, 'Color', SpikeCol(uu,(1:3)))
    LegendPlot{length(SpikeSU)+1-uu} = sprintf('Single Unit %d',uu);
end
set(gca,'YLim',[-1 1]*(NAction+1), 'YTick',[-flip(1:length(SpikeSU))/(length(SpikeSU)+1) 1:NAction], 'YTickLabel', LegendPlot)

orient(Fig,'landscape')
Fig.PaperPositionMode = 'auto';
set(Fig,'Units', 'centimeters', 'Position', get(0, 'screensize'));
set(Fig,'PaperOrientation','landscape');
set(Fig,'PaperUnits','normalized');
set(Fig,'PaperPosition', [0 0 1 1]);

print(Fig,fullfile(Output_dir,sprintf('%s_%s_NeuroActuogram.pdf', Date(3:end), LoggerID)),'-dpdf')
end

               
    


