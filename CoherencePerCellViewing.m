
%% Listing datacells
%Filename = '59834_20190611_SSS_1-97.mat';
% Filename = '59834_20190610_SSS_1-130.mat';
%Path = '/Users/elie/Documents/LMCResults/';
Path = '/Volumes/JulieE8T/LMCResults/';
% Path = '/Users/elie/Documents/ManipBats/LMC/ResultsFiles/';
% Path = '/Users/elie/Google Drive/BatmanData/';

AllFiles = dir(fullfile(Path,'59834*.mat'));
Files2run = zeros(length(AllFiles),1);
for ff=1:length(AllFiles)
    fprintf(1, 'Cell %d/%d\n', ff, length(AllFiles))
    Cell = load(fullfile(AllFiles(ff).folder, AllFiles(ff).name), 'What');
    if (length(strfind(AllFiles(ff).name, '_'))==3) && isfield(Cell, 'What') % Cells without the what field are non-compiled cells probably because they are just noise
        Files2run(ff) = 1;
    end
    clear Cell
end
CellsPath = AllFiles(logical(Files2run));
NCells = length(CellsPath);
%%
for cc=1:NCells
    
    CellPath = fullfile(CellsPath(cc).folder,CellsPath(cc).name);
    fprintf(1, '*** Cell %s %d/%d %s ***\n', CellPath, cc, NCells, FeatureName)
    Cell = load(CellPath, 'What', 'ExpType', 'Who', 'BioSound', 'Duration','SpikesArrivalTimes_Behav','QualitySSU', 'VocOverlap', 'AudioQuality', 'DelayBefore','DelayAfter', 'AuditoryCoherenceFree', 'MotorCoherenceFree','MotorCoherenceOperant');
    if  ~isfield(Cell.MotorCoherenceOperant, 'Error')
        figure(8);
        clf
        subplot(4,1,1)
        histogram(Cell.MotorCoherenceOperant.Bootstrap.Info, 'FaceColor', [0.6350, 0.0780, 0.1840])%red
        ylabel(sprintf('# bootstrap (total = %d)', length(Cell.MotorCoherenceOperant.Bootstrap.Info)))
        xlabel(sprintf('Information on coherence with sound %s', FeatureName))
        hold on
        VL = vline(Cell.MotorCoherenceOperant.Info, '-g');
        VL.LineWidth = 2;
        title(sprintf('Cell Significance of information with bootstraped permutation test p=%.2f', Cell.MotorCoherenceOperant.Info_p))
        
        subplot(4,1,2)
        histogram(Cell.MotorCoherenceOperant.BootstrapFullTime.Info, 'FaceColor', [0.8500, 0.3250, 0.0980])% orange
        ylabel(sprintf('# bootstrap (total = %d)', length(Cell.MotorCoherenceOperant.Bootstrap.Info)))
        xlabel(sprintf('Information on coherence with sound %s', FeatureName))
        hold on
        VL = vline(Cell.MotorCoherenceOperant.Info, '-g');
        VL.LineWidth = 2;
        title(sprintf('Cell Significance of information with bootstraped Full Time shuffling p=%.2f', Cell.MotorCoherenceOperant.Info_pFullTime))
        
        subplot(4,1,3)
        histogram(Cell.MotorCoherenceOperant.BootstrapTime.Info, 'FaceColor', [0.9290, 0.6940, 0.1250]) % yellow
        ylabel(sprintf('# bootstrap (total = %d)', length(Cell.MotorCoherenceOperant.Bootstrap.Info)))
        xlabel(sprintf('Information on coherence with sound %s', FeatureName))
        hold on
        VL = vline(Cell.MotorCoherenceOperant.Info, '-g');
        VL.LineWidth = 2;
        title(sprintf('Cell Significance of information with bootstraped Time shuffling within stims p=%.2f', Cell.MotorCoherenceOperant.Info_pTime))
        
        subplot(4,1,4)
        histogram(Cell.MotorCoherenceOperant.Bootstrap.Info, 'FaceColor', [0.6350, 0.0780, 0.1840])%red
        hold on
        histogram(Cell.MotorCoherenceOperant.BootstrapFullTime.Info, 'FaceColor', [0.8500, 0.3250, 0.0980])% orange
        hold on
        histogram(Cell.MotorCoherenceOperant.BootstrapTime.Info, 'FaceColor', [0.9290, 0.6940, 0.1250]) % yellow
        hold on
        VL = vline(Cell.MotorCoherenceOperant.Info, '-g');
        VL.LineWidth = 2;
        ylabel(sprintf('# bootstrap (total = %d)', length(Cell.MotorCoherenceOperant.Bootstrap.Info)))
        xlabel(sprintf('Information on coherence with sound %s', FeatureName))
        title(sprintf('Cell%d/%d: %s ', cc, NCells,CellsPath(cc).name))
        drawnow
        pause()
    end
end


%%
figure()
subplot(1,3,1)
scatter(MotorCoherenceOperantAll.Info, MotorCoherenceOperantAll.Info_p,20, [0.6350, 0.0780, 0.1840], 'filled')
xlabel('Information (bits)')
ylabel('p-value permutation test')
subplot(1,3,2)
scatter(MotorCoherenceOperantAll.Info, MotorCoherenceOperantAll.Info_pFullTime,20, [0.8500, 0.3250, 0.0980], 'filled')
xlabel('Information (bits)')
ylabel('p-value Time Shuffling Y')
subplot(1,3,3)
scatter(MotorCoherenceOperantAll.Info, MotorCoherenceOperantAll.Info_pTime,20, [0.9290, 0.6940, 0.1250], 'filled')
xlabel('Information (bits)')
ylabel('p-value Time Shuffling Y within Stims')