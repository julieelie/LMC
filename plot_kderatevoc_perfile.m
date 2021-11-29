function plot_kderatevoc_perfile(InputDataFile, OutputPath, Delay )
%%
if nargin<3
    Delay=[3000 200];% in ms
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
% Get the SS quality
NeuralInputID{3} = DataFile(strfind(DataFile, '_SS')+3);

% Get the subject ID
SubjectID = DataFile(1:5);

% Input
FileNameBase = sprintf('%s_%s_SS%s_%s-%s', SubjectID, Date,NeuralInputID{3},NeuralInputID{1},NeuralInputID{2});
FullDataSetFile = fullfile(OutputPath, sprintf('%s.mat', FileNameBase));
load(FullDataSetFile, 'KDE_onset','KDE_offset');


% %% KDE Self vs others all vocalizations from all sessions
% if isfield(KDE_onset, 'OthersVocAll') && isfield(KDE_onset,'SelfVocAll')
%     Fig1 = kdeplot1v1(KDE_onset.SelfVocAll,KDE_onset.OthersVocAll,KDE_offset.SelfVocAll,KDE_offset.OthersVocAll,[1 0 0],[0 0 1],'Self','Others');
%     suplabel('All vocalizations all sessions','t');
%     orient(Fig1,'landscape')
%     print(Fig1,fullfile(OutputPath,sprintf('%s_KDEVoc_AllSession_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')    
% end

%% KDE Trill vs Ba in Operant vs Free, from self
ColOp = [0.8500, 0.3250, 0.0980];
ColFr = [0.6350, 0.0780, 0.1840];
if isfield(KDE_onset, 'SelfTrFr') && isfield(KDE_onset,'SelfBaFr') && isfield(KDE_onset, 'SelfTrOp') && isfield(KDE_onset,'SelfBaOp')
    Fig1 = kdeplot2v2(KDE_onset.SelfTrOp,KDE_onset.SelfBaOp,KDE_onset.SelfTrFr,KDE_onset.SelfBaFr,ColOp,ColFr,'Trill Operant','non-Trill Operant', 'Trill Free', 'non-Trill Free');
    suplabel('Self vocalizations both sessions','t');
    orient(Fig1,'landscape')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfVoc_OpFre_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
end

%% KDE Trill vs Ba in Free from self vs Others
ColOthers = [0, 0.4470, 0.7410];
ColSelf = [0.8500, 0.3250, 0.0980];
if isfield(KDE_onset, 'SelfTrFr') && isfield(KDE_onset,'SelfBaFr') && isfield(KDE_onset, 'OthersTrFr') && isfield(KDE_onset,'OthersBaFr')
    Fig1 = kdeplot2v2(KDE_onset.SelfTrFr,KDE_onset.SelfBaFr,KDE_onset.OthersTrFr,KDE_onset.OthersBaFr,ColSelf,ColOthers,'Trill Self','non-Trill Self', 'Trill Others', 'non-Trill Others');
    suplabel('Vocalizations Free session','t');
    orient(Fig1,'landscape')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfOthersVoc_Fre_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
end

% %% KDE Trill vs Ba from self from all sessions
% BaCol = [0/255 191/255 255/255];
% TrCol = [1, 0, 0];
% if isfield(KDE_onset, 'SelfTrAll') && isfield(KDE_onset,'SelfBaAll')
%     Fig1 = kdeplot1v1(KDE_onset.SelfTrAll,KDE_onset.SelfBaAll,KDE_offset.SelfTrAll,KDE_offset.SelfBaAll,TrCol,BaCol,'Trill','Ba');
%     suplabel('Self vocalizations all sessions','t');
%     orient(Fig1,'landscape')
%     print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfVoc_AllSession_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
% end

%% KDE Trill vs Ba from self from operant sessions
% BaCol = [0/255 191/255 255/255];
% TrCol = [1, 0, 0];
% if isfield(KDE_onset, 'SelfTrOp') && isfield(KDE_onset,'SelfBaOp')
%     Fig1 = kdeplot1v1(KDE_onset.SelfTrOp,KDE_onset.SelfBaOp,KDE_offset.SelfTrOp,KDE_offset.SelfBaOp,TrCol,BaCol,'Trill','Ba');
%     suplabel('Self vocalizations operant session','t');
%     orient(Fig1,'landscape')
%     print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfVoc_OpSession_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
% end

%% KDE Self vocalizations between operant and free
ColOp = [0.8500, 0.3250, 0.0980];
ColFr = [0.6350, 0.0780, 0.1840];
if isfield(KDE_onset, 'SelfVocOp') && isfield(KDE_onset,'SelfVocFr')
    Fig1 = kdeplot1v1_on(KDE_onset.SelfVocOp,KDE_onset.SelfVocFr,ColOp,ColFr,'Operant','Free');
    suplabel('Self vocalizations both sessions','t');
    orient(Fig1,'landscape')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfVoc_OpvsFr_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
end

% %% KDE Self Trills between operant and free
% ColTrFr = [0.75, 0.75, 0];
% TrCol = [0.9290, 0.6940, 0.1250];
% if isfield(KDE_onset, 'SelfTrOp') && isfield(KDE_onset,'SelfTrFr')
%     Fig1 = kdeplot1v1(KDE_onset.SelfTrOp,KDE_onset.SelfTrFr,KDE_offset.SelfTrOp,KDE_offset.SelfTrFr,TrCol,ColTrFr,'Operant','Free session');
%     suplabel('Self Trills all sessions','t');
%     orient(Fig1,'landscape')
%     print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfTr_OpvsFr_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
% end

% %% KDE Self Ba between operant and free
% ColBaFr = [0.4660, 0.6740, 0.1880];
% ColBaOp = [0.4940, 0.1840, 0.5560];
% if isfield(KDE_onset, 'SelfBaOp') && isfield(KDE_onset,'SelfBaFr')
%     Fig1 = kdeplot1v1(KDE_onset.SelfBaOp,KDE_onset.SelfBaFr,KDE_offset.SelfBaOp,KDE_offset.SelfBaFr,ColBaOp,ColBaFr,'Operant','Free session');
%     suplabel('Self Barks all sessions','t');
%     orient(Fig1,'landscape')
%     print(Fig1,fullfile(OutputPath,sprintf('%s_KDESelfBa_OpvsFr_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
% end

%% KDE Others vocalizations between operant and free
% ColOp = [0, 0.4470, 0.7410];
% ColFr = [0.3010, 0.7450, 0.9330];
% if isfield(KDE_onset, 'OthersVocOp') && isfield(KDE_onset,'OthersVocFr')
%     Fig1 = kdeplot1v1(KDE_onset.OthersVocOp,KDE_onset.OthersVocFr,KDE_offset.OthersVocOp,KDE_offset.OthersVocFr,ColOp,ColFr,'Operant','Free session');
%     suplabel('Others vocalizations all sessions','t');
%     orient(Fig1,'landscape')
%     print(Fig1,fullfile(OutputPath,sprintf('%s_KDEOthersVoc_OpvsFr_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
% end

%% KDE all vocalizations Self vs others Operant session
% ColOthersOp = [0, 0.4470, 0.7410];
% ColSelfOp = [0.8500, 0.3250, 0.0980];
% if isfield(KDE_onset, 'OthersVocOp') && isfield(KDE_onset,'SelfVocOp')
%     Fig1 = kdeplot1v1(KDE_onset.OthersVocOp,KDE_onset.SelfVocOp,KDE_offset.OthersVocOp,KDE_offset.SelfVocOp,ColOthersOp,ColSelfOp,'Others','Self');
%     suplabel('All vocalizations Operant session','t');
%     orient(Fig1,'landscape')
%     print(Fig1,fullfile(OutputPath,sprintf('%s_KDEVocOp_SelfvsOthers_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
% end

%% if no comparison possible, just plot the kdeof all vocalizations self in operant session
ColSelfOp = [0.8500, 0.3250, 0.0980];
if isfield(KDE_onset,'SelfVocOp') && ~exist('Fig1','var')
    Fig1 = kdeplot(KDE_onset.SelfVocOp,KDE_offset.SelfVocOp,ColSelfOp);
    suplabel('All vocalizations Operant session','t');
    orient(Fig1,'landscape')
    print(Fig1,fullfile(OutputPath,sprintf('%s_KDEVocOp_Self_%d_%d.pdf', FileNameBase, Delay(1),Delay(2))),'-dpdf','-fillpage')
end

%% Close all figures
close all

%% INTERNAL FUNCTION
    function [FigHand] = kdeplot1v1(Dat1,Dat2,Dat3,Dat4,Col1,Col2,Legend1,Legend2)
        FigHand=figure();
        subplot(1,2,1)
        plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
        hold on
        plot(Dat2(2,:),Dat2(1,:),'-','Color',Col2,'LineWidth',2)
        legend(Legend1,Legend2,'AutoUpdate','off')
        hold on
        shadedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat2(2,:),Dat2(1,:),Dat2(3:4,:),{'-','Color',Col2, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim([-Delay(1) Delay(2)])
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        title('Vocalization Onset')
        
        subplot(1,2,2)
        plot(Dat3(2,:),Dat3(1,:),'-','Color',Col1,'LineWidth',2)
        hold on
        plot(Dat4(2,:),Dat4(1,:),'-','Color',Col2,'LineWidth',2)
        legend(Legend1,Legend2,'AutoUpdate','off')
        hold on
        shadedErrorBar(Dat3(2,:),Dat3(1,:),Dat3(3:4,:),{'-','Color',Col1, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat4(2,:),Dat4(1,:),Dat4(3:4,:),{'-','Color',Col2, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim([-Delay(1) Delay(2)])
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        title('Vocalization Offset')
    end

    function [FigHand] = kdeplot1v1_on(Dat1,Dat2,Col1,Col2,Legend1,Legend2)
        FigHand=figure();
        subplot(1,2,1)
        plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
        hold on
        plot(Dat2(2,:),Dat2(1,:),'-','Color',Col2,'LineWidth',2)
        legend(Legend1,Legend2,'AutoUpdate','off')
        hold on
        shadedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat2(2,:),Dat2(1,:),Dat2(3:4,:),{'-','Color',Col2, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim([-Delay(1) Delay(2)])
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        title('Vocalization Onset')
        
        subplot(1,2,2)
        Xt1 = logical((Dat1(2,:)>=-Delay(1)/10) .* (Dat1(2,:)<=2*Delay(1)/10));
        plot(Dat1(2,Xt1),Dat1(1,Xt1),'-','Color',Col1,'LineWidth',2)
        hold on
        Xt2 = logical((Dat2(2,:)>=-Delay(1)/10) .* (Dat2(2,:)<=2*Delay(1)/10));
        plot(Dat2(2,Xt2),Dat2(1,Xt2),'-','Color',Col2,'LineWidth',2)
        legend(Legend1,Legend2,'AutoUpdate','off')
        hold on
        shadedErrorBar(Dat1(2,Xt1),Dat1(1,Xt1),Dat1(3:4,Xt1),{'-','Color',Col1, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat2(2,Xt2),Dat2(1,Xt2),Dat2(3:4,Xt2),{'-','Color',Col2, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim([-Delay(1)/10 2*Delay(2)/10])
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        title('Vocalization Onset')
    end

    function [FigHand] = kdeplot2v2(Dat1,Dat2,Dat3,Dat4,Col1,Col2,Legend1,Legend2, Legend3, Legend4)
        FigHand=figure();
        subplot(1,2,1)
        plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
        hold on
        plot(Dat2(2,:),Dat2(1,:),'--','Color',Col1,'LineWidth',2)
        hold on
        plot(Dat3(2,:),Dat3(1,:),'-','Color',Col2,'LineWidth',2)
        hold on
        plot(Dat4(2,:),Dat4(1,:),'--','Color',Col2,'LineWidth',2)
        legend(Legend1,Legend2,Legend3, Legend4,'AutoUpdate','off')
        hold on
        shadedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat2(2,:),Dat2(1,:),Dat2(3:4,:),{'--','Color',Col1, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat3(2,:),Dat3(1,:),Dat3(3:4,:),{'-','Color',Col2, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat4(2,:),Dat4(1,:),Dat4(3:4,:),{'--','Color',Col2, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim([-Delay(1) Delay(2)])
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        title('Vocalization Onset')
        
        subplot(1,2,2)
        Xt1 = logical((Dat1(2,:)>=-Delay(1)/10) .* (Dat1(2,:)<=2*Delay(1)/10));
        plot(Dat1(2,Xt1),Dat1(1,Xt1),'-','Color',Col1,'LineWidth',2)
        hold on
        Xt2 = logical((Dat2(2,:)>=-Delay(1)/10) .* (Dat2(2,:)<=2*Delay(1)/10));
        plot(Dat2(2,Xt2),Dat2(1,Xt2),'--','Color',Col1,'LineWidth',2)
        hold on
        Xt3 = logical((Dat3(2,:)>=-Delay(1)/10) .* (Dat3(2,:)<=2*Delay(1)/10));
        plot(Dat3(2,Xt3),Dat3(1,Xt3),'-','Color',Col2,'LineWidth',2)
        hold on
        Xt4 = logical((Dat4(2,:)>=-Delay(1)/10) .* (Dat4(2,:)<=2*Delay(1)/10));
        plot(Dat4(2,Xt4),Dat4(1,Xt4),'--','Color',Col2,'LineWidth',2)
        legend(Legend1,Legend2,Legend3, Legend4,'AutoUpdate','off')
        hold on
        shadedErrorBar(Dat1(2,Xt1),Dat1(1,Xt1),Dat1(3:4,Xt1),{'-','Color',Col1, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat2(2,Xt2),Dat2(1,Xt2),Dat2(3:4,Xt2),{'--','Color',Col1, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat3(2,Xt3),Dat3(1,Xt3),Dat3(3:4,Xt3),{'-','Color',Col2, 'LineWidth',2})
        hold on
        shadedErrorBar(Dat4(2,Xt4),Dat4(1,Xt4),Dat4(3:4,Xt4),{'--','Color',Col2, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim([-Delay(1)/10 2*Delay(2)/10])
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        title('Vocalization Onset')
    end
    

    function [FigHand] = kdeplot(Dat1,Dat3,Col1)
        FigHand=figure();
        subplot(1,2,1)
        plot(Dat1(2,:),Dat1(1,:),'-','Color',Col1,'LineWidth',2)
        hold on
        shadedErrorBar(Dat1(2,:),Dat1(1,:),Dat1(3:4,:),{'-','Color',Col1, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim([-Delay(1) Delay(2)])
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        title('Vocalization Onset')
        
        subplot(1,2,2)
        plot(Dat3(2,:),Dat3(1,:),'-','Color',Col1,'LineWidth',2)
        hold on
        shadedErrorBar(Dat3(2,:),Dat3(1,:),Dat3(3:4,:),{'-','Color',Col1, 'LineWidth',2})
        hold on
        VL = vline(0,':k');
        VL.LineWidth = 2;
        hold off
        xlim([-Delay(1) Delay(2)])
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        title('Vocalization Offset')
    end
end
