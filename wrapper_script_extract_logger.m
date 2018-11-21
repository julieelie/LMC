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
Days = {'20180627', '20180628', '20180702', '20180703', '20180706', '20180710', '20180711', '20180712', '20180724', '20180725', '20180726'};
Server_logger_path = 'Z:\users\JulieE\GiMo_65430_71300\Loggers';
Server_audio_path = 'Z:\users\JulieE\GiMo_65430_71300\Audio';
for dd=1:length(Days)
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

%%  