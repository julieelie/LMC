%% Gather files of operant and free session that have been manually curated
BasePath = '/Volumes/JulieE8T/LMC_CoEd/';
[List2AudioPath, SessionType] = gather_audio_datapath(BasePath);
Nsets = length(List2AudioPath);

%% Investigate some statistics (nb of calls...)
Ncalls = 0;
NTrills = 0;
BatIDs = cell(Nsets,1);
Mean_SpecMean_piezo = cell(Nsets,1);
Mean_Fund_piezo = cell(Nsets,1);
Mean_Amp_piezo = cell(Nsets,1);
Mean_Saliency_piezo = cell(Nsets,1);
SessionID = cell(Nsets,1);
for Seti=1:Nsets
    load(List2AudioPath{Seti}, 'BioSoundCalls', 'BioSoundFilenames');
    NVoc = size(BioSoundCalls,1);
    BatIDs{Seti} = cell(NVoc,1);
    Mean_SpecMean_piezo{Seti} = nan(NVoc,1);
    Mean_Fund_piezo{Seti} = nan(NVoc,1);
    Mean_Amp_piezo{Seti} = nan(NVoc,1);
    Mean_Saliency_piezo{Seti} = nan(NVoc,1);
    if strcmp(SessionType{Seti}, 'O')
        SessionID{Seti} = ones(NVoc,1);
    else
        SessionID{Seti} = zeros(NVoc,1);
    end
    for vv=1:NVoc
        Ind_AL = strfind(BioSoundFilenames{vv,1},'_AL');
        Ind_Bat = strfind(BioSoundFilenames{vv,1}, '_Bat');
        BatIDs{Seti}{vv} = BioSoundFilenames{vv,1}((Ind_Bat + 5):(Ind_AL-1));
        Mean_SpecMean_piezo{Seti}(vv) = nanmean(BioSoundCalls{vv,2}.SpectralMean);
        Mean_Fund_piezo{Seti}(vv) = nanmean(BioSoundCalls{vv,2}.f0);
        Mean_Amp_piezo{Seti}(vv) = nanmean(BioSoundCalls{vv,2}.SoundAmp);
        Mean_Saliency_piezo{Seti}(vv) = nanmean(BioSoundCalls{vv,2}.sal);
    end 
    Ncalls = NVoc + Ncalls;
    NTrills = sum(contains(BioSoundCalls(:,1).type, 'Tr')) + NTrills;
end 

BatIDs = [BatIDs{:}];
Mean_SpecMean_piezo = [Mean_SpecMean_piezo{:}];
Mean_Fund_piezo = [Mean_Fund_piezo{:}];
Mean_Amp_piezo = [Mean_Amp_piezo{:}];
Mean_Saliency_piezo = [Mean_Saliency_piezo{:}];
SessionID = [SessionID{:}];


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
%ExpFolders = dir(fullfile(BasePath,'LMC*'));
ExpFolders = BasePath;
NF = 0; % counter for single files
for ee=1:length(ExpFolders)
    fprintf(1, '\n  -> Looking into  %s...\n ', fullfile(ExpFolders(ee).folder,ExpFolders(ee).name))
    DateFolders = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name, 'logger','20*'));
    for dd=1:length(DateFolders)
        fprintf(1, '   %s\n', DateFolders(dd).name);
        AudioDataFiles = dir(fullfile(DateFolders(dd).folder, DateFolders(dd).name,'*_VocEXtractData_*'));
        for ll = 1:length(AudioDataFiles)
            NF = NF +1;
            List2ParamPath{NF} = fullfile(AudioDataFiles(ll).folder, AudioDataFiles(ll).name);
            % finding session type
            ParamFile = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name,'audio',DateFolders(dd).name, sprintf('%s_%s*param.txt',ExpFolders(ee).folder(5:8) , AudioDataFiles(ll).name(1:11))));
            if contains(ParamFile, 'VocTrigger')
                SessionType{NF} = 'O';
            elseif contains(ParamFile, 'RecOnly')
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
