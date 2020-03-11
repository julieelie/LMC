function wrapper_who_calls_check(Path2ParamFile, Logger_dir)
addpath(genpath('/Users/elie/Documents/CODE/LMC'))
addpath(genpath('/Users/elie/Documents/CODE/LoggerDataProcessing'))
addpath(genpath('/Users/elie/Documents/CODE/SoundAnalysisBats'))
close all
% Get the recording data
[AudioDataPath, DataFile ,~]=fileparts(Path2ParamFile);
Date = DataFile(6:11);

if nargin<2
    % Set the path to logger data
    Logger_dir = fullfile(AudioDataPath(1:(strfind(AudioDataPath, 'audio')-1)), 'logger',['20' Date]);
end



% Checking what we have in terms of vocalization localization/extraction
ExpStartTime = DataFile(13:16);


%% Identify who is calling
fprintf(' CHECKING RESULT OF IDENTIFY WHO IS CALLING\n')
who_calls_check(AudioDataPath, Logger_dir,Date, ExpStartTime,200,1);


end