function neuralData_compile_perfile(InputDataFile, OutputPath, MaxDur)

% INPUT:
%       InputDataFile: path to a single unit mat file.

%       OutputPath: string that indicate where the data should be saved. If
%       not specified or left empty, the data will be saved in the folder
%       of the inputfile
if nargin<3
    MaxDur = 500; % Duration in ms for estimating the rate
end

% Get the paths
[Path2Data, DataFile]=fileparts(InputDataFile);
PathParts = regexp(Path2Data, '/', 'split');
Loggers_dir = ['/' fullfile(PathParts{1:end-2})];
if nargin<2 || isempty(OutputPath)
    OutputPath = Loggers_dir;
end

% Get the date of the recording
Idx_ = strfind(DataFile, '_');
Date = DataFile((Idx_(1)+1) : (Idx_(2)-1));

% Get the tetrode ID
NeuralInputID{1} = DataFile(strfind(DataFile, 'TT')+2);
% Get the SS ID
NeuralInputID{2} = DataFile((Idx_(end)+1):end);

% Get the subject ID
SubjectID = DataFile(1:5);

%% Get the logger ID for each bat

%% Load the other behavior actions
% Initialize output matrix
DataDir = dir(fullfile(OutputPath, sprintf('%s_%s_*_SSU%s-%s.mat', SubjectID, Date,NeuralInputID{1},NeuralInputID{2})));
for ff=1:length(DataDir)
    Neuro = load(fullfile(DataDir(ff).folder, DataDir(ff).name));
    FieldNames = fieldnames(Neuro);
    for fi=1:length(FieldNames)
        if strcmp(FieldNames{fi}, 'Voc_NeuroSSU')
            Idx_2 = strfind(DataDir(ff).name,'_');
            ExpStartTime = DataDir(ff).name((Idx_2(2) +1):(Idx_2(3)-1));
            DataFile = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData_*.mat', Date(3:end), ExpStartTime)));
            load(fullfile(DataFile.folder, DataFile.name), 'IndVocStartRaw_merged', 'IndVocStopRaw_merged', 'BatID','LoggerName');
            load(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractData.mat', Date(3:end), ExpStartTime)), 'FS','Piezo_wave','Raw_wave');
            
            
        end
    end
end
Num_Slots(bb) = sum(floor(Dur_local/MaxDur)) + sum(mod(Dur_local, MaxDur)>0);