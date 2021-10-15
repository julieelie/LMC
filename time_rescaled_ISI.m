function [ISIvecnew, ISIvecold, Rescaled_spikes]=time_rescaled_ISI(Spiketrains, RateVecs)
% this function gets time rescaled ISIs according to time varying rate RateVecs and original ISIs ISI = InterSpike intervals
% SpikeTrains is a cell array of spike patterns sample at 1000Hz with one
% cell per stimulus/vocalization/event
% RateVecs is the estimated time varying rate (FS = 1000Hz) for each of
% these stimuli/vocalizations/events
DebugFig=0;
NEvents=size(Spiketrains, 1);
ISIvecnew=cell(NEvents,1);
ISIvecold=cell(NEvents,1);
Burstvec=[.45 .3 .2];  
Rescaled_spikes = cell(size(Spiketrains));
if NEvents==1 && ~iscell(RateVecs)
    RateVecs = {RateVecs};
end
for ee=1:NEvents
    Meanrate=mean(RateVecs{ee});
    if isnan(Meanrate) || Meanrate==0
        fprintf(1,'no spikes in mean for Event %d\n', ee);
        Rescaled_spikes{ee}=Spiketrains{ee};
        continue
    end
    Rescaled_spikes{ee}=zeros(size(Spiketrains{ee}));
    Time_indices=find(Spiketrains{ee});
    Extratimes=[];
    %the following takes care of bins with more than one spike. Extra
    %spikes are spread in the same ms bin. if 2 spikes then the first is at
    %t and second at t+0.45, if 3 spikes then they are placed at t, t+0.3,
    %t+0.6; if 4 spikes then they are placed at t, t+0.2, t+0.4 and t+0.6
    if any(Spiketrains{ee}>1)
        for Nspikes=2:4
            Burstindex=find(Spiketrains{ee}==Nspikes);
            for ns=1:Nspikes-1
                Extratimes=[Extratimes Burstindex+Burstvec(Nspikes-1)*ns];
            end
        end
    elseif any(Spiketrains{ee}>4)
        warning('The code does not handle more than 4 spikes per bin')
        keyboard
    end
    % To remove the time avrying rate from the underlying spike noise
    % time between spikes is strechted (high rate) or condensed (low rate)
    Newtime=cumsum(RateVecs{ee})/Meanrate;
    Ntime_indices=Newtime(Time_indices); % Here are the rate corrected spike arrival times
    % these lines take care of the burst cases (more than one spike per
    % bin) by finding the corrected time arrival of these spikes using
    % cubic interpolation in Newtime
    if ~isempty(Extratimes)
        Ntimeindicesextras=spline(1:length(RateVecs{ee}), Newtime, Extratimes);
        Ntime_indices=sort([Ntime_indices Ntimeindicesextras]);
    end
    
    if DebugFig
        figure(18);
        clf
        subplot(2,1,1)
        hist(Time_indices,length(Spiketrains{ee}))
        title(sprintf('NSpike = %d',sum(Spiketrains{ee})))
        xlim([0 length(Spiketrains{ee})])
        ylabel('Original spike pattern')
        h=findobj(gca, 'Type', 'patch');
        h.FaceColor = 'k';
        h.EdgeColor = 'k';
        hold on
        yyaxis right
        plot(RateVecs{ee}*10^3,'b-', 'LineWidth',2)
        hold on
        plot(Newtime, 'b--', 'LineWidth',2)
        ylabel('Original rate Hz and New time')
        hold off
        subplot(2,1,2)
        hist(Ntime_indices,length(Spiketrains{ee}))
        title(sprintf('NSpike = %d',sum(Spiketrains{ee})))
        xlim([0 length(Spiketrains{ee})])
        ylabel('Rescaled spike pattern')
        h=findobj(gca, 'Type', 'patch');
        h.FaceColor = 'k';
        h.EdgeColor = 'k';
        hold on
        yyaxis right
        plot(RateVecs{ee}*10^3,'b-', 'LineWidth',2)
        hold on
        plot(Newtime, 'b--', 'LineWidth',2)
        ylabel('Original rate mHz and New time')
        hold off
        pause()
    end
    
    ISIvecnew{ee}=diff(Ntime_indices);
    ISIvecold{ee}=diff(sort([Time_indices Extratimes]));
    if ~isempty(Ntime_indices)
        %this doesn't make a one if there is a zero but does not deal with
        % duplicates
        Rescaled_spikes{ee}(round(Ntime_indices(round(Ntime_indices)~=0)))=1;
        if length(unique(round(Ntime_indices))) ~= length(Ntime_indices)
            warning('%d Spikes are arriving in the same bin and are lost in when constructing the rescaled spike trains', length(Ntime_indices)-length(unique(round(Ntime_indices))))
%             keyboard
        end
    end
    
end