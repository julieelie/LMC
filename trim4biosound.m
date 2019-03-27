function []=trim4biosound(Loggers_dir, Date, ExpStartTime, MergeThresh, Manual)
% Loggers_dir = '/Users/elie/Documents/ManipBats/LMC/190110_59882_11689_HoHa/20190130_TestCalls';
Loggers_dir = '/Users/elie/Documents/ManipBats/JulieBatsDrive/180711/loggers';
%Date = num2str(190130);
Date = num2str(180711);
%ExpStartTime = num2str(1007);
ExpStartTime = num2str(1059);
%MergeThresh=200; %ms
BandPassFiltPiezo = [100 10000]; % bandpass filter of the piezo

%% Load the data
load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)))
% load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date, ExpStartTime)), 'Piezo_wave', 'Piezo_FS',  'Raw_wave','FS','BandPassFilter', 'AudioLogs','VocFilename');
% load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_%d.mat', Date, ExpStartTime, MergeThresh)), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'IndVocStartPiezo_merged', 'IndVocStopPiezo_merged');

%% Creating the output directories that will contain the wavfiles
OutputDir = fullfile(Loggers_dir, '4Biosound');
OutputDirMic = fullfile(OutputDir, 'Mic');
OutputDirPiezo = fullfile(OutputDir, 'Piezo');
mkdir(OutputDir)
mkdir(OutputDirMic)
mkdir(OutputDirPiezo)

%% Create correctly cut files for biosound
Nsound = length(Raw_wave);
LoggerNames = fieldnames(Piezo_wave);
for ss=1:Nsound
    fprintf(1,'Treating sound %d/%d\n', ss, Nsound)
    % Check that only and exactly one bat was vocalizing
    VocBatInd=find(~cellfun('isempty', IndVocStartRaw_merged{ss}));
    if length(VocBatInd)>1 || isempty(VocBatInd)
        continue
    end
    % get the name of that file
    [~,File]=fileparts(VocFilename{ss});
    % cut the mic and piezo files along the proposed indices
    Ncuts = length(IndVocStartRaw_merged{ss}{VocBatInd});
    NcutsP = length(IndVocStartPiezo_merged{ss}{VocBatInd});
    if NcutsP~=Ncuts
        error('There is not the same number of cuts for sound %d/%d in the piezo and Raw records\n',ss,Nsound)
    end
    for cc=1:Ncuts
        fprintf(1,'    ->Mic\n')
        Y = Raw_wave{ss}(IndVocStartRaw_merged{ss}{VocBatInd}(cc):IndVocStopRaw_merged{ss}{VocBatInd}(cc));
        Filename_local = fullfile(OutputDirMic, sprintf('%s_mic_%d.wav',File,cc));
        Y=Y-mean(Y);% center the sound
        audiowrite(Filename_local, Y, FS);
        fprintf(1,'    ->Piezo\n')
        Y = Piezo_wave.(LoggerNames{VocBatInd}){ss}(IndVocStartPiezo_merged{ss}{VocBatInd}(cc):IndVocStopPiezo_merged{ss}{VocBatInd}(cc));
        [z,p,k] = butter(6,BandPassFiltPiezo/(Piezo_FS.(LoggerNames{VocBatInd})(ss)/2),'bandpass');
        sos_low = zp2sos(z,p,k);
        Y = Y-mean(Y); % center the output of the piezo
        Y = Y/max(abs(Y)); % scale it between -1 and 1
        YLow = (filtfilt(sos_low,1,Y)); % low-pass filter the voltage trace
        Filename_local = fullfile(OutputDirPiezo, sprintf('%s_piezo_%d_%s.wav',File,cc,LoggerNames{VocBatInd}));
        audiowrite(Filename_local, YLow/max(abs(YLow)), round(Piezo_FS.(LoggerNames{VocBatInd})(ss)));
    end
        
end

end