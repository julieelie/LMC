%% Gather files of operant and free session that have been manually curated
BasePath = '/Volumes/JulieE8T/';
[List2AudioPath, SessionType] = gather_audio_datapath(BasePath);
Nsets = length(List2AudioPath);

%% Gather data
WarningID = 'MATLAB:Python:UnsupportedLoad';
warning('off', WarningID)
Ncalls = 0;
NTrills = 0;
BatIDs = cell(1,Nsets);
Mean_SpecMean_piezo = cell(1,Nsets);
Mean_Fund_piezo = cell(1,Nsets);
Mean_Amp_piezo = cell(1,Nsets);
Mean_Saliency_piezo = cell(1,Nsets);
Duration_ms = cell(1,Nsets);
SessionID = cell(1,Nsets);
Trills01 = cell(1,Nsets);
Filenames = cell(1,Nsets);
for Seti=1:Nsets
    clear BioSoundFilenames
    fprintf(1,'Set %d/%d %s\n', Seti, Nsets, List2AudioPath{Seti})
    load(List2AudioPath{Seti}, 'BioSoundCalls', 'BioSoundFilenames');
    if exist('BioSoundFilenames', 'var')
        VocInd = find(~cellfun('isempty',(BioSoundFilenames(:,1))));
        NVoc = length(VocInd);
        Filenames{Seti} = cell(1,NVoc);
        BatIDs{Seti} = cell(1,NVoc);
        Mean_SpecMean_piezo{Seti} = nan(1,NVoc);
        Mean_Fund_piezo{Seti} = nan(1,NVoc);
        Mean_Amp_piezo{Seti} = nan(1,NVoc);
        Mean_Saliency_piezo{Seti} = nan(1,NVoc);
        Duration_ms{Seti} = nan(1,NVoc);
        Trills01{Seti} = nan(1,NVoc);
        if strcmp(SessionType{Seti}, 'O')
            SessionID{Seti} = ones(1,NVoc);
        else
            SessionID{Seti} = zeros(1,NVoc);
        end
        for ii=1:NVoc
            vv=VocInd(ii);
            Ind_AL = strfind(BioSoundFilenames{vv,1},'_AL');
            Ind_Bat = strfind(BioSoundFilenames{vv,1}, '_Bat');
            BatIDs{Seti}{ii} = BioSoundFilenames{vv,1}((Ind_Bat + 4):(Ind_AL-1));
            Mean_SpecMean_piezo{Seti}(ii) = nanmean(BioSoundCalls{vv,2}.SpectralMean(~isnan(BioSoundCalls{vv,2}.sal)));
            if isnan(Mean_SpecMean_piezo{Seti}(ii))
                keyboard
            end
            Mean_Fund_piezo{Seti}(ii) = nanmean(BioSoundCalls{vv,2}.f0);
            Mean_Amp_piezo{Seti}(ii) = nanmean(BioSoundCalls{vv,2}.amp(~isnan(BioSoundCalls{vv,2}.sal)));
            Mean_Saliency_piezo{Seti}(ii) = BioSoundCalls{vv,2}.meansal;
            Trills01{Seti}(ii) = contains(BioSoundCalls{vv,1}.type, 'Tr');
            Duration_ms{Seti}(ii) = length(BioSoundCalls{vv,2}.amp);
            Filenames{Seti}{ii} = BioSoundFilenames{vv,2};
        end 
        Ncalls = NVoc + Ncalls;
        NTrills = sum(Trills01{Seti}) + NTrills;
        clear BioSoundFilenames
    end
end 
warning('on', WarningID)
BatIDs = [BatIDs{:}]';
Mean_SpecMean_piezo = [Mean_SpecMean_piezo{:}]';
Mean_Fund_piezo = [Mean_Fund_piezo{:}]';
Mean_Amp_piezo = [Mean_Amp_piezo{:}]';
Mean_Saliency_piezo = [Mean_Saliency_piezo{:}]';
SessionID = [SessionID{:}]';
Trills01 = [Trills01{:}]';
Filenames = [Filenames{:}]';
Duration_ms = [Duration_ms{:}]';

%% Investigate some statistics (nb of calls...)

% Number of calls per bat
BatIDx = unique(str2double(BatIDs));
CountByID = nan(length(BatIDx),2);
for bb=1:length(BatIDx)
    CountByID(bb,1) = sum(str2double(BatIDs(logical(SessionID)))==BatIDx(bb));
    CountByID(bb,2) = sum(str2double(BatIDs(~SessionID))==BatIDx(bb));
end
figure();
FIG1 = subplot(1,2,1)
BAR = bar(CountByID);
FIG1.XTickLabel = unique(BatIDs);
title('Total number of calls over 12 free sessions and 33 operant sessions')
xlabel('BatID')
ylabel('# calls')
FIG2 = subplot(1,2,2);
BAR2 = bar(CountByID./([33.*ones(6,1) 12.*ones(6,1)]));
FIG2.XTickLabel = unique(BatIDs);
title('Average number of calls per session (33 OP 12Free)')
xlabel('BatID')
ylabel('# calls')

%% Acoustic landscape of calls in both sessions
figure(2)
scatter(Mean_Amp_piezo, Duration_ms,20,[0 0 0])
xlabel('Amplitude')
ylabel('Duration (ms)')
title(sprintf('n=%d',sum(~isnan(Mean_Amp_piezo))))

Short = Duration_ms<500;

figure(2)
scatter(Mean_Amp_piezo(Short), Duration_ms(Short),20,[0 0 0])
xlabel('Amplitude')
ylabel('Duration (ms)')
title(sprintf('n=%d',sum(~isnan(Mean_Amp_piezo(Short)))))

%%
figure(1)

subplot(2,3,1)
scatter(Mean_Amp_piezo, Mean_Saliency_piezo,10,[0 0 0])
xlabel('Amplitude')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo))))

subplot(2,3,2)
scatter(Mean_Amp_piezo, Mean_SpecMean_piezo,10,[0 0 0])
xlabel('Amplitude')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo))))

subplot(2,3,3)
scatter(Mean_Amp_piezo, Mean_Fund_piezo,10,[0 0 0])
xlabel('Amplitude')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo))))

subplot(2,3,4)
scatter(Duration_ms(Short), Mean_Saliency_piezo(Short),10,[0 0 0])
xlabel('Duration (ms)')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(Short)))))

subplot(2,3,5)
scatter(Duration_ms(Short), Mean_SpecMean_piezo(Short),10,[0 0 0])
xlabel('Duration (ms)')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(Short)))))

subplot(2,3,6)
scatter(Duration_ms(Short), Mean_Fund_piezo(Short),10,[0 0 0])
xlabel('Duration (ms)')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(Short)))))
suplabel('The soundscape of all calls')

%% How are trills dfferent from the rest?
figure(1)

subplot(2,3,1)
scatter(Mean_Amp_piezo, Mean_Saliency_piezo,10,Trills01*[0.4660, 0.6740, 0.1880])
xlabel('Amplitude')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo))))

subplot(2,3,2)
scatter(Mean_Amp_piezo, Mean_SpecMean_piezo,10,Trills01*[0.4660, 0.6740, 0.1880])
xlabel('Amplitude')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo))))

subplot(2,3,3)
scatter(Mean_Amp_piezo, Mean_Fund_piezo,10,Trills01*[0.4660, 0.6740, 0.1880])
xlabel('Amplitude')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo))))

subplot(2,3,4)
scatter(Duration_ms(Short), Mean_Saliency_piezo(Short),10,Trills01(Short)*[0.4660, 0.6740, 0.1880])
xlabel('Duration (ms)')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(Short)))))

subplot(2,3,5)
scatter(Duration_ms(Short), Mean_SpecMean_piezo(Short),10,Trills01(Short)*[0.4660, 0.6740, 0.1880])
xlabel('Duration (ms)')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(Short)))))

subplot(2,3,6)
scatter(Duration_ms(Short), Mean_Fund_piezo(Short),10,Trills01(Short)*[0.4660, 0.6740, 0.1880])
xlabel('Duration (ms)')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(Short)))))
suplabel('Trills in this soundscape')

%% Are sessions covering a different sound landscape?
figure(1)
subplot(2,3,1)
scatter(Mean_Amp_piezo, Mean_Saliency_piezo,10,SessionID*[0, 0.4470, 0.7410] + (~SessionID)*[0.8500, 0.3250, 0.0980])
xlabel('Amplitude')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo))))

subplot(2,3,2)
scatter(Mean_Amp_piezo, Mean_SpecMean_piezo,10,SessionID*[0, 0.4470, 0.7410] + (~SessionID)*[0.8500, 0.3250, 0.0980])
xlabel('Amplitude')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo))))

subplot(2,3,3)
scatter(Mean_Amp_piezo, Mean_Fund_piezo,10,SessionID*[0, 0.4470, 0.7410] + (~SessionID)*[0.8500, 0.3250, 0.0980])
xlabel('Amplitude')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo))))

subplot(2,3,4)
scatter(Duration_ms(Short), Mean_Saliency_piezo(Short),10,SessionID(Short)*[0, 0.4470, 0.7410] + (~SessionID(Short))*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(Short)))))

subplot(2,3,5)
scatter(Duration_ms(Short), Mean_SpecMean_piezo(Short),10,SessionID(Short)*[0, 0.4470, 0.7410] + (~SessionID(Short))*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(Short)))))

subplot(2,3,6)
scatter(Duration_ms(Short), Mean_Fund_piezo(Short),10,SessionID(Short)*[0, 0.4470, 0.7410] + (~SessionID(Short))*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(Short)))))
suplabel('Session Type in this soundscape')

%% How is Cooper covering that sound landscape?
figure(1)
 
subplot(2,3,1)
scatter(Mean_Amp_piezo, Mean_Saliency_piezo,10,(str2double(BatIDs)==59834)*[0.3010, 0.7450, 0.9330])
xlabel('Amplitude')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo))))

subplot(2,3,2)
scatter(Mean_Amp_piezo, Mean_SpecMean_piezo,10,(str2double(BatIDs)==59834)*[0.3010, 0.7450, 0.9330])
xlabel('Amplitude')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo))))

subplot(2,3,3)
scatter(Mean_Amp_piezo, Mean_Fund_piezo,10,(str2double(BatIDs)==59834)*[0.3010, 0.7450, 0.9330])
xlabel('Amplitude')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo))))

subplot(2,3,4)
scatter(Duration_ms(Short), Mean_Saliency_piezo(Short),10,(str2double(BatIDs(Short))==59834)*[0.3010, 0.7450, 0.9330])
xlabel('Duration (ms)')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(Short)))))

subplot(2,3,5)
scatter(Duration_ms(Short), Mean_SpecMean_piezo(Short),10,(str2double(BatIDs(Short))==59834)*[0.3010, 0.7450, 0.9330])
xlabel('Duration (ms)')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(Short)))))

subplot(2,3,6)
scatter(Duration_ms(Short), Mean_Fund_piezo(Short),10,(str2double(BatIDs(Short))==59834)*[0.3010, 0.7450, 0.9330])
xlabel('Duration (ms)')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(Short)))))

suplabel('Calls of Cooper in this soundscape')
%% How is Edwardo covering that sound landscape?
figure(1)

subplot(2,3,1)
scatter(Mean_Amp_piezo, Mean_Saliency_piezo,10,(str2double(BatIDs)==65701)*[0.3010, 0.7450, 0.9330])
xlabel('Amplitude')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo))))

subplot(2,3,2)
scatter(Mean_Amp_piezo, Mean_SpecMean_piezo,10,(str2double(BatIDs)==65701)*[0.3010, 0.7450, 0.9330])
xlabel('Amplitude')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo))))

subplot(2,3,3)
scatter(Mean_Amp_piezo, Mean_Fund_piezo,10,(str2double(BatIDs)==65701)*[0.3010, 0.7450, 0.9330])
xlabel('Amplitude')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo))))

subplot(2,3,4)
scatter(Duration_ms(Short), Mean_Saliency_piezo(Short),10,(str2double(BatIDs(Short))==65701)*[0.3010, 0.7450, 0.9330])
xlabel('Duration (ms)')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(Short)))))

subplot(2,3,5)
scatter(Duration_ms(Short), Mean_SpecMean_piezo(Short),10,(str2double(BatIDs(Short))==65701)*[0.3010, 0.7450, 0.9330])
xlabel('Duration (ms)')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(Short)))))

subplot(2,3,6)
scatter(Duration_ms(Short), Mean_Fund_piezo(Short),10,(str2double(BatIDs(Short))==65701)*[0.3010, 0.7450, 0.9330])
xlabel('Duration (ms)')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(Short)))))

suplabel('Calls of Edwardo in this soundscape')

%% Are coopers calls different between operant and free session?
SessionID_Cooper = SessionID(str2double(BatIDs)==59834);
SessionID_cooperShort = SessionID(logical(Short.*str2double(BatIDs)==59834));
figure(1)

subplot(2,3,1)
scatter(Mean_Amp_piezo(str2double(BatIDs)==59834), Mean_Saliency_piezo(str2double(BatIDs)==59834),10,SessionID_Cooper*[0, 0.4470, 0.7410] + (~SessionID_Cooper)*[0.8500, 0.3250, 0.0980])
xlabel('Amplitude')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(str2double(BatIDs)==59834)))))

subplot(2,3,2)
scatter(Mean_Amp_piezo(str2double(BatIDs)==59834), Mean_SpecMean_piezo(str2double(BatIDs)==59834),10,SessionID_Cooper*[0, 0.4470, 0.7410] + (~SessionID_Cooper)*[0.8500, 0.3250, 0.0980])
xlabel('Amplitude')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(str2double(BatIDs)==59834)))))

subplot(2,3,3)
scatter(Mean_Amp_piezo(str2double(BatIDs)==59834), Mean_Fund_piezo(str2double(BatIDs)==59834),10,SessionID_Cooper*[0, 0.4470, 0.7410] + (~SessionID_Cooper)*[0.8500, 0.3250, 0.0980])
xlabel('Amplitude')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(str2double(BatIDs)==59834)))))

subplot(2,3,4)
scatter(Duration_ms(logical(Short.*str2double(BatIDs)==59834)), Mean_Saliency_piezo(logical(Short.*str2double(BatIDs)==59834)),10,SessionID_cooperShort*[0, 0.4470, 0.7410] + (~SessionID_cooperShort)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(logical(Short.*str2double(BatIDs)==59834))))))

subplot(2,3,5)
scatter(Duration_ms(logical(Short.*str2double(BatIDs)==59834)), Mean_SpecMean_piezo(logical(Short.*str2double(BatIDs)==59834)),10,SessionID_cooperShort*[0, 0.4470, 0.7410] + (~SessionID_cooperShort)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(logical(Short.*str2double(BatIDs)==59834))))))

subplot(2,3,6)
scatter(Duration_ms(logical(Short.*str2double(BatIDs)==59834)), Mean_Fund_piezo(logical(Short.*str2double(BatIDs)==59834)),10,SessionID_cooperShort*[0, 0.4470, 0.7410] + (~SessionID_cooperShort)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(logical(Short.*str2double(BatIDs)==59834))))))


suplabel('Calls of Cooper different between Free and Operant?')

%% Are Trill calls different between operant and free session?
SessionID_Trills = SessionID(logical(Trills01));
SessionID_Trills_Short = SessionID(logical(Trills01 .* Short));
figure(1)

subplot(2,3,1)
scatter(Mean_Amp_piezo(logical(Trills01)), Mean_Saliency_piezo(logical(Trills01)),10,SessionID_Trills*[0, 0.4470, 0.7410] + (~SessionID_Trills)*[0.8500, 0.3250, 0.0980])
xlabel('Amplitude')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(logical(Trills01))))))
hold on
MeanAmp_OP = mean(Mean_Amp_piezo(logical(Trills01 .*SessionID)));
MeanAmp_FS = mean(Mean_Amp_piezo(logical(Trills01 .*~SessionID)));
MeanSal_OP = mean(Mean_Saliency_piezo(logical(Trills01 .*SessionID)));
MeanSal_FS = mean(Mean_Saliency_piezo(logical(Trills01 .*~SessionID)));
plot(MeanAmp_OP, MeanSal_OP,'o','MarkerSize',15,'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerEdgeColor',[0, 0, 0])
hold on
plot(MeanAmp_FS, MeanSal_FS,'o','MarkerSize',15,'MarkerFaceColor',[0.8500, 0.3250, 0.0980], 'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
hold off

subplot(2,3,2)
scatter(Mean_Amp_piezo(logical(Trills01)), Mean_SpecMean_piezo(logical(Trills01)),10,SessionID_Trills*[0, 0.4470, 0.7410] + (~SessionID_Trills)*[0.8500, 0.3250, 0.0980])
MeanMeanSpecMean_OP = mean(Mean_SpecMean_piezo(logical(Trills01 .*SessionID)));
MeanMeanSpecMean_FS = mean(Mean_SpecMean_piezo(logical(Trills01 .*~SessionID)));
hold on
plot(MeanAmp_OP, MeanMeanSpecMean_OP,'o','MarkerSize',15,'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerEdgeColor',[0, 0, 0])
hold on
plot(MeanAmp_FS, MeanMeanSpecMean_FS,'o','MarkerSize',15,'MarkerFaceColor',[0.8500, 0.3250, 0.0980], 'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
hold off
xlabel('Amplitude')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(logical(Trills01))))))

subplot(2,3,3)
scatter(Mean_Amp_piezo(logical(Trills01)), Mean_Fund_piezo(logical(Trills01)),10,SessionID_Trills*[0, 0.4470, 0.7410] + (~SessionID_Trills)*[0.8500, 0.3250, 0.0980])
MeanMeanFund_OP = mean(Mean_Fund_piezo(logical(Trills01 .*SessionID)));
MeanMeanFund_FS = nanmean(Mean_Fund_piezo(logical(Trills01 .*~SessionID)));
hold on
plot(MeanAmp_OP, MeanMeanFund_OP,'o','MarkerSize',15,'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerEdgeColor',[0, 0, 0])
hold on
plot(MeanAmp_FS, MeanMeanFund_FS,'o','MarkerSize',15,'MarkerFaceColor',[0.8500, 0.3250, 0.0980], 'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
hold off
xlabel('Amplitude')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(logical(Trills01))))))

subplot(2,3,4)
scatter(Duration_ms(logical(Short.*Trills01)), Mean_Saliency_piezo(logical(Short.*Trills01)),10,SessionID_Trills_Short*[0, 0.4470, 0.7410] + (~SessionID_Trills_Short)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(logical(Short.*Trills01))))))

subplot(2,3,5)
scatter(Duration_ms(logical(Short.*Trills01)), Mean_SpecMean_piezo(logical(Short.*Trills01)),10,SessionID_Trills_Short*[0, 0.4470, 0.7410] + (~SessionID_Trills_Short)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(logical(Short.*Trills01))))))

subplot(2,3,6)
scatter(Duration_ms(logical(Short.*Trills01)), Mean_Fund_piezo(logical(Short.*Trills01)),10,SessionID_Trills_Short*[0, 0.4470, 0.7410] + (~SessionID_Trills_Short)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(logical(Short.*Trills01))))))


suplabel('Trill calls from Cooper different between Free and Operant?')

%% Longer Louder?
figure()
scatter(Duration_ms(logical(Short.*Trills01)), Mean_Amp_piezo(logical(Short.*Trills01)),10,SessionID_Trills_Short*[0, 0.4470, 0.7410] + (~SessionID_Trills_Short)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Amplitude')
title(sprintf('n=%d',sum(~isnan(Mean_Amp_piezo(logical(Short.*Trills01))))))

%% Are non-Trill calls different between operant and free session?
SessionID_Barks_Cooper = SessionID(logical(~Trills01.*(str2double(BatIDs)==59834)));
SessionID_Barks_Short_Cooper = SessionID(logical((~Trills01) .* Short.*(str2double(BatIDs)==59834)));
figure(1)

subplot(2,3,1)
scatter(Mean_Amp_piezo(logical(~Trills01.*(str2double(BatIDs)==59834))), Mean_Saliency_piezo(logical(~Trills01.*(str2double(BatIDs)==59834))),10,SessionID_Barks_Cooper*[0, 0.4470, 0.7410] + (~SessionID_Barks_Cooper)*[0.8500, 0.3250, 0.0980])
xlabel('Amplitude')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(logical(~Trills01.*(str2double(BatIDs)==59834)))))))
hold on
MeanAmp_OP = nanmean(Mean_Amp_piezo(logical((~Trills01.*(str2double(BatIDs)==59834)) .*SessionID)));
MeanAmp_FS = nanmean(Mean_Amp_piezo(logical((~Trills01.*(str2double(BatIDs)==59834)) .*~SessionID)));
MeanSal_OP = nanmean(Mean_Saliency_piezo(logical((~Trills01.*(str2double(BatIDs)==59834)) .*SessionID)));
MeanSal_FS = nanmean(Mean_Saliency_piezo(logical((~Trills01.*(str2double(BatIDs)==59834)) .*~SessionID)));
plot(MeanAmp_OP, MeanSal_OP,'o','MarkerSize',15,'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerEdgeColor',[0, 0, 0])
hold on
plot(MeanAmp_FS, MeanSal_FS,'o','MarkerSize',15,'MarkerFaceColor',[0.8500, 0.3250, 0.0980], 'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
hold off

subplot(2,3,2)
scatter(Mean_Amp_piezo(logical(~Trills01.*(str2double(BatIDs)==59834))), Mean_SpecMean_piezo(logical(~Trills01.*(str2double(BatIDs)==59834))),10,SessionID_Barks_Cooper*[0, 0.4470, 0.7410] + (~SessionID_Barks_Cooper)*[0.8500, 0.3250, 0.0980])
MeanMeanSpecMean_OP = nanmean(Mean_SpecMean_piezo(logical((~Trills01.*(str2double(BatIDs)==59834)) .*SessionID)));
MeanMeanSpecMean_FS = nanmean(Mean_SpecMean_piezo(logical((~Trills01.*(str2double(BatIDs)==59834)) .*~SessionID)));
hold on
plot(MeanAmp_OP, MeanMeanSpecMean_OP,'o','MarkerSize',15,'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerEdgeColor',[0, 0, 0])
hold on
plot(MeanAmp_FS, MeanMeanSpecMean_FS,'o','MarkerSize',15,'MarkerFaceColor',[0.8500, 0.3250, 0.0980], 'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
hold off
xlabel('Amplitude')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(logical(~Trills01.*(str2double(BatIDs)==59834)))))))

subplot(2,3,3)
scatter(Mean_Amp_piezo(logical(~Trills01.*(str2double(BatIDs)==59834))), Mean_Fund_piezo(logical(~Trills01.*(str2double(BatIDs)==59834))),10,SessionID_Barks_Cooper*[0, 0.4470, 0.7410] + (~SessionID_Barks_Cooper)*[0.8500, 0.3250, 0.0980])
MeanMeanFund_OP = nanmean(Mean_Fund_piezo(logical((~Trills01.*(str2double(BatIDs)==59834)) .*SessionID)));
MeanMeanFund_FS = nanmean(Mean_Fund_piezo(logical((~Trills01.*(str2double(BatIDs)==59834)) .*~SessionID)));
hold on
plot(MeanAmp_OP, MeanMeanFund_OP,'o','MarkerSize',15,'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerEdgeColor',[0, 0, 0])
hold on
plot(MeanAmp_FS, MeanMeanFund_FS,'o','MarkerSize',15,'MarkerFaceColor',[0.8500, 0.3250, 0.0980], 'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
hold off
xlabel('Amplitude')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(logical(~Trills01.*(str2double(BatIDs)==59834)))))))

subplot(2,3,4)
scatter(Duration_ms(logical(Short.*~Trills01.*(str2double(BatIDs)==59834))), Mean_Saliency_piezo(logical(Short.*~Trills01.*(str2double(BatIDs)==59834))),10,SessionID_Barks_Short_Cooper*[0, 0.4470, 0.7410] + (~SessionID_Barks_Short_Cooper)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Picth saliency')
title(sprintf('n=%d',sum(~isnan(Mean_Saliency_piezo(logical(Short.*(~Trills01).*(str2double(BatIDs)==59834)))))))

subplot(2,3,5)
scatter(Duration_ms(logical(Short.*(~Trills01).*(str2double(BatIDs)==59834))), Mean_SpecMean_piezo(logical(Short.*(str2double(BatIDs)==59834).*(~Trills01))),10,SessionID_Barks_Short_Cooper*[0, 0.4470, 0.7410] + (~SessionID_Barks_Short_Cooper)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Spectral Mean')
title(sprintf('n=%d',sum(~isnan(Mean_SpecMean_piezo(logical(Short.*(~Trills01).*(str2double(BatIDs)==59834)))))))

subplot(2,3,6)
scatter(Duration_ms(logical(Short.*(str2double(BatIDs)==59834).*~Trills01)), Mean_Fund_piezo(logical(Short.*(~Trills01).*(str2double(BatIDs)==59834))),10,SessionID_Barks_Short_Cooper*[0, 0.4470, 0.7410] + (~SessionID_Barks_Short_Cooper)*[0.8500, 0.3250, 0.0980])
xlabel('Duration (ms)')
ylabel('Fundamental')
title(sprintf('n=%d',sum(~isnan(Mean_Fund_piezo(logical(Short.*(~Trills01).*(str2double(BatIDs)==59834)))))))


suplabel('Bark calls from Cooper different between Free and Operant?')


%% Construct Trills vector
Trills = zeros(Ncalls,1);
for vv=1:Ncalls
    Trills(vv) = strcmp(BioSoundCalls{vv,1}.type, 'Tr');
end

%% Plot the Periodicity
AmpPeriodPThresh = 0.01;
AmpPeriod = nan(Ncalls,2);
Fund = nan(Ncalls,1);
ColorCode = nan(Ncalls,3);
FigC=figure()
subplot(2,2,1)
for vv=1:Ncalls
    if ~isempty(BioSoundCalls{vv,1}.AmpPeriodF)
        AmpPeriod(vv,1) = BioSoundCalls{vv,1}.AmpPeriodF;
        AmpPeriod(vv,2) = BioSoundCalls{vv,1}.AmpPeriodP;
        try
            Fund(vv) = BioSoundCalls{vv,2}.fund;
        catch
            Fund(vv) = NaN;
        end
    end
    hold on
    if Trills(vv)
        ColorCode(vv,:) = [1 0 0];
        plot(BioSoundCalls{vv,1}.AmpPeriodF, BioSoundCalls{vv,1}.AmpPeriodP, 'ro','MarkerFaceColor', 'r','MarkerSize',10,'MarkerEdgeColor','r');
    else
        ColorCode(vv,:) = [0 0 0];
        plot(BioSoundCalls{vv,1}.AmpPeriodF, BioSoundCalls{vv,1}.AmpPeriodP, 'ro','MarkerFaceColor', 'k','MarkerSize',10,'MarkerEdgeColor','k');
    end
end
ylabel('Amplitude Periodicity Index Mic')
xlabel('Amplitude Periodicity Frequency Mic (Hz)')

subplot(2,2,2)
for vv=1:Ncalls
    if BioSoundCalls{vv,1}.AmpPeriodP>AmpPeriodPThresh
        hold on
        if Trills(vv)
            plot(BioSoundCalls{vv,1}.AmpPeriodF, BioSoundCalls{vv,1}.AmpPeriodP, 'ro','MarkerFaceColor', 'r','MarkerSize',10,'MarkerEdgeColor','r');
        else
            plot(BioSoundCalls{vv,1}.AmpPeriodF, BioSoundCalls{vv,1}.AmpPeriodP, 'ro','MarkerFaceColor', 'k','MarkerSize',10,'MarkerEdgeColor','k');
        end
    end
end
ylabel('Amplitude Periodicity Index Mic')
xlabel('Amplitude Periodicity Frequency Mic (Hz)')
subplot(2,2,3)
for vv=1:Ncalls
    if ~isempty(BioSoundCalls{vv,2}.AmpPeriodF)
        AmpPeriod(vv,1) = BioSoundCalls{vv,2}.AmpPeriodF;
        AmpPeriod(vv,2) = BioSoundCalls{vv,2}.AmpPeriodP;
    end
    hold on
    if Trills(vv)
        plot(BioSoundCalls{vv,2}.AmpPeriodF, BioSoundCalls{vv,2}.AmpPeriodP, 'ro','MarkerFaceColor', 'r','MarkerSize',10,'MarkerEdgeColor','r');
    else
        plot(BioSoundCalls{vv,2}.AmpPeriodF, BioSoundCalls{vv,2}.AmpPeriodP, 'ro','MarkerFaceColor', 'k','MarkerSize',10,'MarkerEdgeColor','k');
    end
end
ylabel('Amplitude Periodicity Index Piezo')
xlabel('Amplitude Periodicity Frequency Piezo (Hz)')

subplot(2,2,4)
for vv=1:Ncalls
    if BioSoundCalls{vv,2}.AmpPeriodP>AmpPeriodPThresh
        hold on
        if Trills(vv)
            plot(BioSoundCalls{vv,2}.AmpPeriodF, BioSoundCalls{vv,2}.AmpPeriodP, 'ro','MarkerFaceColor', 'r','MarkerSize',10,'MarkerEdgeColor','r');
        else
            plot(BioSoundCalls{vv,2}.AmpPeriodF, BioSoundCalls{vv,2}.AmpPeriodP, 'ro','MarkerFaceColor', 'k','MarkerSize',10,'MarkerEdgeColor','k');
        end
    end
end
ylabel('Amplitude Periodicity Index Piezo')
xlabel('Amplitude Periodicity Frequency Piezo (Hz)')
print(FigC,fullfile(Loggers_dir,'SoundAnalysis','Periodicity.pdf'),'-dpdf','-fillpage')

%% Plot periodicity and fundamental
FigD = figure()
scatter3(AmpPeriod(:,1),AmpPeriod(:,2),Fund,20,ColorCode,'filled')
xlabel('Amplitude Periodicity Frequency Mic (Hz)')
ylabel('Amplitude Periodicity Index Mic')
zlabel('Fundamental Frequency (Hz)')
print(FigD,fullfile(Loggers_dir,'SoundAnalysis','PeriodicityFund.pdf'),'-dpdf','-fillpage')

%% Plot Saliency and Spectral mean
FigB=figure()
for vv=1:Ncalls
    hold on
    if Trills(vv)
        plot(BioSoundCalls{vv,1}.meansal, nanmean(BioSoundCalls{vv,1}.SpectralMean), 'ro','MarkerFaceColor', 'r','MarkerSize',10,'MarkerEdgeColor','r');
    else
        plot(BioSoundCalls{vv,1}.meansal, nanmean(BioSoundCalls{vv,1}.SpectralMean), 'ro','MarkerFaceColor', 'k','MarkerSize',10,'MarkerEdgeColor','k');
    end
    
end
ylabel('Spectral Mean (Hz)')
xlabel('Pitch Saliency')
print(FigB,fullfile(Loggers_dir,'SoundAnalysis','SalSpecMean.pdf'),'-dpdf','-fillpage')

%% Temporal parameters  (std time, entropy time)
FigA=figure()
subplot(1,2,1)
for vv=1:Ncalls
    hold on
    if Trills(vv)
        plot(BioSoundCalls{vv,1}.stdtime, nanmean(BioSoundCalls{vv,1}.entropytime), 'ro','MarkerFaceColor', 'r','MarkerSize',10,'MarkerEdgeColor','r');
    else
        plot(BioSoundCalls{vv,1}.stdtime, nanmean(BioSoundCalls{vv,1}.entropytime), 'ro','MarkerFaceColor', 'k','MarkerSize',10,'MarkerEdgeColor','k');
    end
    
end
ylabel('Duration (std amplitude enveloppe)')
xlabel('Temporal struct (entropy amplitude enveloppe)')

subplot(1,2,2)
for vv=1:Ncalls
    hold on
    if Trills(vv)
        plot(BioSoundCalls{vv,1}.stdtime, nanmean(BioSoundCalls{vv,1}.entropytime), 'ro','MarkerFaceColor', 'r','MarkerSize',10,'MarkerEdgeColor','r');
    else
        plot(BioSoundCalls{vv,1}.stdtime, nanmean(BioSoundCalls{vv,1}.entropytime), 'ro','MarkerFaceColor', 'k','MarkerSize',10,'MarkerEdgeColor','k');
    end
    
end
xlim([0 0.15])
ylabel('Duration (std amplitude enveloppe)')
xlabel('Temporal struct (entropy amplitude enveloppe)')
print(FigA,fullfile(Loggers_dir,'SoundAnalysis','StdTimeEntropyTime.pdf'),'-dpdf','-fillpage')

%% Plot the dynamic of all calls
Fig=figure()
for vv=1:Ncalls
    hold on
    Legend = (vv==1);
    plotCallDynamic(BioSoundCalls{vv,1}, BioSoundCalls{vv,2}, Legend)
    
end
%title ('All Calls')
print(Fig,fullfile(Loggers_dir,'SoundAnalysis','Dynamic_AllCalls.pdf'),'-dpdf','-fillpage')

Fig=figure()
IndTr = find(Trills);
for vv=1:length(IndTr)
    hold on
    Legend = (vv==1);
    plotCallDynamic(BioSoundCalls{IndTr(vv),1}, BioSoundCalls{IndTr(vv),2}, Legend)
    
end
%title ('Trills')
print(Fig,fullfile(Loggers_dir,'SoundAnalysis','Dynamic_Trills.pdf'),'-dpdf','-fillpage')

Fig=figure()
IndBa = find(~Trills);
for vv=1:length(IndBa)
    hold on
    Legend = (vv==1);
    plotCallDynamic(BioSoundCalls{IndBa(vv),1}, BioSoundCalls{IndBa(vv),2}, Legend)
    
end
%title ('Barks')
print(Fig,fullfile(Loggers_dir,'SoundAnalysis','Dynamic_Barks.pdf'),'-dpdf','-fillpage')

%% INTERNAL FUNCTION

function [List2ParamPath,SessionType] = gather_audio_datapath(BasePath)
fprintf(1,'*** Gathering paths to audio reconly data ***')
List2ParamPath = cell(10^3,1); % initialize the list to 1000
SessionType = cell(10^3,1); % initialize the list to 1000
ExpFolders = dir(fullfile(BasePath,'LMC*'));
%ExpFolders = BasePath;
NF = 0; % counter for single files
for ee=1:length(ExpFolders)
    fprintf(1, '\n  -> Looking into  %s...\n ', fullfile(ExpFolders(ee).folder,ExpFolders(ee).name))
    DateFolders = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name, 'logger','20*'));
    for dd=1:length(DateFolders)
        fprintf(1, '   %s\n', DateFolders(dd).name);
        AudioDataFiles = dir(fullfile(DateFolders(dd).folder, DateFolders(dd).name,'*_VocExtractDat*_*'));
        for ll = 1:length(AudioDataFiles)
            NF = NF +1;
            List2ParamPath{NF} = fullfile(AudioDataFiles(ll).folder, AudioDataFiles(ll).name);
            % finding session type
            ParamFile = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name,'audio',DateFolders(dd).name, sprintf('%s_%s*param.txt',ExpFolders(ee).name(5:8) , AudioDataFiles(ll).name(1:11))));
            if contains(ParamFile.name, 'VocTrigger')
                SessionType{NF} = 'O';
            elseif contains(ParamFile.name, 'RecOnly')
                SessionType{NF} = 'F';
            end
        end
    end
end
List2ParamPath = List2ParamPath(1:NF);
SessionType = SessionType(1:NF);
fprintf(1, '\n %d files of operant sessions and %d files of free recording sessions have been retrieved\n', sum(contains(SessionType, 'O')), sum(contains(SessionType, 'F')));
end

function plotCallDynamic(BiosoundRaw, BiosoundPiezo,Legend)
Span = 9;% Span is an unevennumber. smooth has a default span of 5 points = 5ms However end points are unchanged...
HalfSpan = (Span-1)/2;
YLIM = [0 0.4];
% Plot the pitch saliency vs amplitude on microphone
subplot(4,1,1)
Saliency = mysmooth(double(BiosoundRaw.sal), Span);
TimeSound = double(BiosoundRaw.to)*1000;
TimeSound = TimeSound./max(TimeSound);
cmap = colormap('jet');
ncolors = length(cmap);
nx = length(Saliency);

for ii=HalfSpan:nx-HalfSpan
    segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
    plot([Saliency(ii), Saliency(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
    hold on;
end
set(gca,'XLim',[0 1]);
set(gca,'YLim',YLIM);
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


for ii=HalfSpan:nx-HalfSpan
    segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
    plot([FormantDisp(ii), FormantDisp(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
    hold on;
end
set(gca,'XLim',[10 150])
set(gca,'YLim',YLIM)
if Legend
    xlabel('1/Formant disp (vocal tract length (mm))')
    ylabel('Amplitude')
end

% Plot the amplitude (Mic data) vs fundamental (Piezo
% Data)
subplot(4,1,3)
SoundFund = mysmooth(double(BiosoundPiezo.f0), Span);
if ~isempty(SoundFund)
    for ii=HalfSpan:nx-HalfSpan
        segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
        plot([SoundFund(ii), SoundFund(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
        hold on;
    end
    set(gca,'YLim',YLIM)
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
    for ii=HalfSpan:nx-HalfSpan
        segcolor = cmap(fix((TimeSound(ii)+TimeSound(ii+1))*ncolors./3)+1,:);
        plot([SoundSpecMean(ii), SoundSpecMean(ii+1)], [BiosoundRaw.amp(ii), BiosoundRaw.amp(ii+1)], "Color",segcolor, "LineWidth",2);
        hold on;
    end
    set(gca,'YLim',YLIM)
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
