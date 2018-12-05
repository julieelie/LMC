%% Path to code
addpath(genpath('C:\Users\Julie\Documents\GitHub\GeneralCode'))
addpath(genpath('C:\Users\Julie\Documents\GitHub\LMC'))
addpath(genpath('C:\Users\Julie\Documents\GitHub\LoggerDataProcessing'))
addpath(genpath('C:\Users\Julie\Documents\GitHub\SoundAnalysisBats'))
Days = {'20180627', '20180628', '20180702', '20180703', '20180706', '20180710', '20180711', '20180712', '20180724', '20180725', '20180726'};
Server_logger_path = 'Z:\users\JulieE\GiMo_65430_71300\Loggers';
Server_audio_path = 'Z:\users\JulieE\GiMo_65430_71300\Audio';
VocDir = 'Z:\users\JulieE\GiMo_65430_71300\Audio\ManualCuts';
NeuroLoggerID = 'Logger16';
Path2NlxCode = 'C:\Users\Julie\Documents\GitHub\MatlabImportExport_v6.0.0';

%% Convert the spike sorted ntt files to mat files

for dd=1:length(Days)
    fprintf(' CONVERTING NTT FILES TO MAT FILES \n')
    fprintf('*********** %s *************\n', Days{dd})
    Loggers_dir = fullfile(Server_logger_path, Days{dd});
    InputPath = fullfile(Loggers_dir,NeuroLoggerID,'extracted_data');
    ntt2mat(InputPath,Path2NlxCode);
end

%% Extract the audio logger data of the days we are interested in
script_extract_loggers('06272018',{'65430' '71300' '71354' '71335' '71332'} , [5 6 3 7 9])
script_extract_loggers('06282018',{'65430' '71300' '71354' '71213' '71137'} , [5 6 3 7 9])
script_extract_loggers('07022018',{'65430' '65430' '71300' '71300' '59899'} , [5 7 9 6 3])
script_extract_loggers('07032018',{'65430' '65430' '71300' '71300' '71354'} , [5 7 9 6 3])
script_extract_loggers('07062018',{'65430' '71300' '60014'} , [5 6 3])
script_extract_loggers('07102018',{'65430' '71300' '59899' '71132' '71213'} , [5 6 10 9 7])
script_extract_loggers('07112018',{'65430' '71300' '59899' '71132' '71335'} , [5 6 10 9 7])
script_extract_loggers('07122018',{'65430' '71300' '71173' '71137' '71213'} , [5 6 10 9 7])
script_extract_loggers('07242018',{'65430' '71300' '59899' '71132' '71335'} , [5 6 10 9 7])
script_extract_loggers('07252018',{'65430' '71300' '60014' '71132' '71213'} , [5 6 10 9 7])
script_extract_loggers('07262018',{'65430' '71300' '71173' '71137' '71335'} , [5 6 10 9 7])

%% Allign TTL pulses for all these days

for dd=10:length(Days)
    fprintf('*********** %s *************\n', Days{dd})
    Loggers_dir = fullfile(Server_logger_path, Days{dd});
    Audio_dir = fullfile(Server_audio_path, Days{dd});
    % find the time stamp of each experiment
    StampFiles = dir(fullfile(Audio_dir, '*RecOnly_param.txt'));
    if length(StampFiles)>1
        fprintf('Several Recording Only tests were done on that day, please choose the one you want to look at:\n')
        for ff=1:length(StampFiles)
            fprintf('%d: %s\n',ff, StampFiles(ff).name);
        end
        Indff = input('Your choice:\n');
    else
        Indff = 1;
    end
    ExpStartTime = StampFiles(Indff).name(13:16);
    align_soundmexAudio_2_logger(Audio_dir, Loggers_dir, ExpStartTime,'Method','rise');
end

%% Now for each day and each manualy extracted vocalization, find the index of the corresponding first sample in the environmental microphone recordings

for dd=1:length(Days)
    fprintf(' LOCALIZING VOCALIZATIONS ON AMBIENT MIC \n')
    fprintf('*********** %s *************\n', Days{dd})
    Audio_dir = fullfile(Server_audio_path, Days{dd});
    % find the time stamp of each experiment
    StampFiles = dir(fullfile(Audio_dir, '*RecOnly_param.txt'));
    if length(StampFiles)>1
        fprintf('Several Recording Only tests were done on that day, please choose the one you want to look at:\n')
        for ff=1:length(StampFiles)
            fprintf('%d: %s\n',ff, StampFiles(ff).name);
        end
        Indff = input('Your choice:\n');
    else
        Indff = 1;
    end
    ExpStartTime = StampFiles(Indff).name(13:16);
    voc_localize(VocDir, Audio_dir, Days{dd}(3:end), ExpStartTime);
end

%% Identify the same vocalizations on the piezos and save sound extracts, onset and offset times

for dd=1:length(Days)
    fprintf(' LOCALIZING VOCALIZATIONS ON PIEZO RECORDINGS\n')
    fprintf('*********** %s *************\n', Days{dd})
    Audio_dir = fullfile(Server_audio_path, Days{dd});
    Loggers_dir = fullfile(Server_logger_path, Days{dd});
    % find the time stamp of each experiment
    StampFiles = dir(fullfile(Audio_dir, '*RecOnly_param.txt'));
    if length(StampFiles)>1
        fprintf('Several Recording Only tests were done on that day, please choose the one you want to look at:\n')
        for ff=1:length(StampFiles)
            fprintf('%d: %s\n',ff, StampFiles(ff).name);
        end
        Indff = input('Your choice:\n');
    else
        Indff = 1;
    end
    ExpStartTime = StampFiles(Indff).name(13:16);
    get_logger_data_voc(Audio_dir, Loggers_dir,Days{dd}(3:end), ExpStartTime);
end
%% Identify who is calling
for dd=1:length(Days)
    fprintf(' IDENTIFY WHO IS CALLING\n')
    fprintf('*********** %s *************\n', Days{dd})
    Audio_dir = fullfile(Server_audio_path, Days{dd});
    Loggers_dir = fullfile(Server_logger_path, Days{dd});
    % find the time stamp of each experiment
    StampFiles = dir(fullfile(Audio_dir, '*RecOnly_param.txt'));
    if length(StampFiles)>1
        fprintf('Several Recording Only tests were done on that day, please choose the one you want to look at:\n')
        for ff=1:length(StampFiles)
            fprintf('%d: %s\n',ff, StampFiles(ff).name);
        end
        Indff = input('Your choice:\n');
    else
        Indff = 1;
    end
    ExpStartTime = StampFiles(Indff).name(13:16);
    who_calls(Loggers_dir,Days{dd}(3:end), ExpStartTime,200,0);
end

%% Extract the neural data corresponding to the vocalizations
for dd=1:length(Days)
    fprintf(' EXTRACTING NEURAL DATA CORRESPONDING TO VOCALIZATIONS \n')
    fprintf('*********** %s *************\n', Days{dd})
    Audio_dir = fullfile(Server_audio_path, Days{dd});
    Loggers_dir = fullfile(Server_logger_path, Days{dd});
    % find the time stamp of each experiment
    StampFiles = dir(fullfile(Audio_dir, '*RecOnly_param.txt'));
    if length(StampFiles)>1
        fprintf('Several Recording Only tests were done on that day, please choose the one you want to look at:\n')
        for ff=1:length(StampFiles)
            fprintf('%d: %s\n',ff, StampFiles(ff).name);
        end
        Indff = input('Your choice:\n');
    else
        Indff = 1;
    end
    ExpStartTime = StampFiles(Indff).name(13:16);
    cut_neuralData_voc(Loggers_dir,Days{dd}(3:end), ExpStartTime,[0 1 1 1],200);
end


%% Plot PSTH of the bats hearing or producing a vocalization
addpath(genpath('C:\Users\Julie\Documents\GitHub\GeneralCode\'))
AudioLoggerID = {'Logger5';'Logger5' ; 'Logger7' ;'Logger7';'Logger5';'Logger5';'Logger5';'Logger5';'Logger5';'Logger5';'Logger5'};
NeuroLoggerID = 'Logger16';
Flags=[1 1];
for dd=1:length(Days)
    fprintf(' PSTH of NEURAL DATA CORRESPONDING TO VOCALIZATIONS \n')
    fprintf('*********** %s *************\n', Days{dd})
    Audio_dir = fullfile(Server_audio_path, Days{dd});
    Loggers_dir = fullfile(Server_logger_path, Days{dd});
    % find the time stamp of each experiment
    StampFiles = dir(fullfile(Audio_dir, '*RecOnly_param.txt'));
    if length(StampFiles)>1
        fprintf('Several Recording Only tests were done on that day, please choose the one you want to look at:\n')
        for ff=1:length(StampFiles)
            fprintf('%d: %s\n',ff, StampFiles(ff).name);
        end
        Indff = input('Your choice:\n');
    else
        Indff = 1;
    end
    ExpStartTime = StampFiles(Indff).name(13:16);
    plot_psth_voc(Loggers_dir, Days{dd}(3:end), ExpStartTime, AudioLoggerID{dd}, NeuroLoggerID, Flags, 200)
    close all
%     pause()
end

%% Plot PSTH of the bats doing other actions!
addpath(genpath('C:\Users\Julie\Documents\GitHub\GeneralCode\'))
AudioLoggerID = {'Logger5';'Logger5' ; 'Logger7' ;'Logger7';'Logger5';'Logger5';'Logger5';'Logger5';'Logger5';'Logger5';'Logger5'};
NeuroLoggerID = 'Logger16';
Flags=[1 1];
for dd=1:length(Days)
    fprintf(' PSTH of NEURAL DATA CORRESPONDING TO ALL BEHAVIORS \n')
    fprintf('*********** %s *************\n', Days{dd})
    Audio_dir = fullfile(Server_audio_path, Days{dd});
    Loggers_dir = fullfile(Server_logger_path, Days{dd});
    % find the time stamp of each experiment
    StampFiles = dir(fullfile(Audio_dir, '*RecOnly_param.txt'));
    if length(StampFiles)>1
        fprintf('Several Recording Only tests were done on that day, please choose the one you want to look at:\n')
        for ff=1:length(StampFiles)
            fprintf('%d: %s\n',ff, StampFiles(ff).name);
        end
        Indff = input('Your choice:\n');
    else
        Indff = 1;
    end
    ExpStartTime = StampFiles(Indff).name(13:16);
    ExtData_dir = fullfile(Loggers_dir, NeuroLoggerID, 'extracted_data');
    get_logger_data_event(ExtData_dir, Days{dd}, 200, NeuroLoggerID)
    close all
%     pause()
end

