function plot_kderatevoc_perfile(InputDataFile, OutputPath)

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
FileNameBase = sprintf('%s_%s_SSU%s-%s', SubjectID, Date,NeuralInputID{1},NeuralInputID{2});
FullDataSetFile = fullfile(OutputPath, sprintf('%s.mat', FileNameBase));
load(FullDataSetFile, 'KDE');


%% KDE Self vs others all vocalizations from all sessions
if isfield(KDE, 'OthersVocAll') && isfield(KDE,'SelfVocAll')
    Fig1 = kdeplot1v1(KDE.SelfVocAll,KDE.OthersVocAll,[1 0 0],[0 0 1],'Self','Others');
    title('All vocalizations all sessions')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDEVoc_AllSession.pdf', FileNameBase)),'-dpdf','-fillpage')
end

%% KDE Trill vs Ba from self from all sessions
BaCol = [0.9290, 0.6940, 0.1250];
if isfield(KDE, 'SelfTrAll') && isfield(KDE,'SelfBaAll')
    Fig1 = kdeplot1v1(KDE.SelfTrAll,KDE.SelfBaAll,[1 0 0],BaCol,'Trill','Ba');
    title('Self vocalizations all sessions')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfVoc_AllSession.pdf', FileNameBase)),'-dpdf','-fillpage')
end

%% KDE Self vocalizations between operant and free
ColOp = [0.8500, 0.3250, 0.0980];
ColFr = [0.6350, 0.0780, 0.1840];
if isfield(KDE, 'SelfVocOp') && isfield(KDE,'SelfVocFr')
    Fig1 = kdeplot1v1(KDE.SelfVocOp,KDE.SelfVocFr,ColOp,ColFr,'Operant','Free session');
    title('Self vocalizations all sessions')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfVoc_OpvsFr.pdf', FileNameBase)),'-dpdf','-fillpage')
end

%% KDE Self Trills between operant and free
ColTrFr = [0.75, 0.75, 0];
if isfield(KDE, 'SelfTrOp') && isfield(KDE,'SelfTrFr')
    Fig1 = kdeplot1v1(KDE.SelfTrOp,KDE.SelfTrFr,[1 0 0],ColTrFr,'Operant','Free session');
    title('Self Trills all sessions')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfTr_OpvsFr.pdf', FileNameBase)),'-dpdf','-fillpage')
end

%% KDE Self Ba between operant and free
ColBaFr = [0.4660, 0.6740, 0.1880];
ColBaOp = [0.4940, 0.1840, 0.5560];
if isfield(KDE, 'SelfBaOp') && isfield(KDE,'SelfBaFr')
    Fig1 = kdeplot1v1(KDE.SelfBaOp,KDE.SelfBaFr,ColBaOp,ColBaFr,'Operant','Free session');
    title('Self Barks all sessions')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfBa_OpvsFr.pdf', FileNameBase)),'-dpdf','-fillpage')
end

%% KDE Others vocalizations between operant and free
ColOp = [0, 0.4470, 0.7410];
ColFr = [0.3010, 0.7450, 0.9330];
if isfield(KDE, 'OthersVocOp') && isfield(KDE,'OthersVocFr')
    Fig1 = kdeplot1v1(KDE.OthersVocOp,KDE.OthersVocFr,ColOp,ColFr,'Operant','Free session');
    title('Others vocalizations all sessions')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDEOthersVoc_OpvsFr.pdf', FileNameBase)),'-dpdf','-fillpage')
end

%% INTERNAL FUNCTION
    function [FigHand] = kdeplot1v1(Dat1,Dat2,Col1,Col2,Legend1,Legend2)
        FigHand=figure();
        plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
        hold on
        plot(Dat2(2,:),Dat2(1,:),'-','Color',Col2,'LineWidth',2)
        legend(Legend1,Legend2,'AutoUpdate','off')
        hold on
        shaddedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
        hold on
        shaddedErrorBar(Dat2(2,:),Dat2(1,:),Dat2(3:4,:),{'-','Color',Col2, 'LineWidth',2})
        hold off
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
    end
    
end
