% BasePath = '/Volumes/JulieE8T';
BasePath = '/Volumes/server_home/users/JulieE/LMC/';
BaseDataDir = '/Volumes/server_home/users/JulieE/DeafSalineGroup151/';
BaseCodeDir = '/Users/elie/Documents/CODE/GitHub/';
%% Paths to code
addpath(genpath(fullfile(BaseCodeDir,'LMC')))
addpath(genpath(fullfile(BaseCodeDir, 'LoggerDataProcessing')))
addpath(genpath(fullfile(BaseCodeDir,'SoundAnalysisBats')))
addpath(genpath(fullfile(BaseCodeDir,'General')))

%% Get the list of recording sessions
List2RecOnlyPath = gather_reconly_datapath(BasePath);

Path2RunRecOnly = 1:length(List2RecOnlyPath);
% Path2RunRecOnly = find(contains(List2RecOnlyPath, 'CoEd'));
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '201905'))=[];% No neural data
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '190605_1553'))=[]; % Clock jump, no data to extract
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190703'))=[];% No TTL Pulses
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190709'))=[];% No TTL Pulses
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190619'))=[];% Issue of clock drift for logger 49 and 12

% Path2RunRecOnly = find(contains(List2RecOnlyPath, 'HoHa'));
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190131'))=[];% Allignment issue as of now no neural data extracted
Path2RunRecOnly(contains(List2RecOnlyPath(Path2RunRecOnly), '20190202_1400'))=[];% No logger data

%%
% Error = nan(length(Path2RunRecOnly),1);
for pp= 35:length(Path2RunRecOnly)
    Path2ParamFile = List2RecOnlyPath{Path2RunRecOnly(pp)};
    [Folder,File,~] = fileparts(Path2ParamFile);
    Date = File(6:11);
    if str2double(Date)<190116
        Error(pp)=0;
        continue % There's no Free session before that date
    end
    ExpStartTime = File(13:16);
    Logger_dir = fullfile(Folder(1:(strfind(Folder, 'audio')-1)), 'logger');
    % apply BatIDNameCorrect
    fprintf(1, '***********************************************\n* Running BatIDNameCorrect on:\n %s\n Date: %s\n ExpStartTime:%s *\n***********************************************\n', Logger_dir, Date, ExpStartTime)
    Error(pp) = BatIDNameCorrect(Logger_dir, Date, ExpStartTime);
end

%% sum Error=0 so no correction is necessary

%% INTERNAL FUNCTION
function Error = BatIDNameCorrect(Loggers_dir, Date, ExpStartTime)
Error=0; %counter for the number of vocalizations misslabeled
% Find the google drive folder
GGFolder = '/Users/elie/Google Drive/My Drive/';
if ~(exist(GGFolder, 'file')==7)
    GGFolder = '/Users/elie/Google Drive/Mon Drive/';
end
if ~(exist(GGFolder, 'file')==7)
    warning('cannot find GGFolder at %s\n', GGFolder)
    GGFolder = input('Please enter the GGFolder path: ', 's');
    keyboard
end

% Set the path to the recording log
if contains(Loggers_dir, 'LMC')
    Path2RecordingTable = fullfile(GGFolder,'/BatmanData/RecordingLogs/recording_logs.xlsx');
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:P200','basic');
    RowData = find((cell2mat(RecTableData(2:end,1))== str2double(Date))) +1;
elseif contains(Loggers_dir, 'Juvenile')
    Path2RecordingTable = fullfile(GGFolder,'/JuvenileRecordings/JuvenileRecordingsNWAF155_Log.xlsx');
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:Q123','basic');
    RowData = find((cell2mat(RecTableData(4:end,1))== str2double(['20' Date]))) +3;
elseif contains(Loggers_dir, 'DeafSaline')
    Path2RecordingTable = fullfile(GGFolder,'/JuvenileRecordings/DeafRecordingsNWAF155_Log.xlsx');
    [~,~,RecTableData]=xlsread(Path2RecordingTable,1,'A1:AE41','basic');
    RowData = find((cell2mat(RecTableData(2:end,1))== str2double(['20' Date]))) +1;
end
% Find the ID of the recorded bats
DataInfo = RecTableData(RowData,:);
Header = RecTableData(1,:);
BatIDCol = find(contains(Header, 'Bat'));
All_loggers_dir = dir(fullfile(Loggers_dir,['20' Date], '*ogger*'));
DirFlags = [All_loggers_dir.isdir];
% Extract only those that are directories.
All_loggers_dir = All_loggers_dir(DirFlags);
LoggerName = cell(length(All_loggers_dir),1);
BatID = cell(length(All_loggers_dir),1);
for ll=1:length(All_loggers_dir)
    Ind = strfind(All_loggers_dir(ll).name, 'r');
    Logger_num = str2double(All_loggers_dir(ll).name((Ind+1):end));
    NLogCol = find(contains(Header, 'NL'));% Columns of the neural loggers
    ALogCol = find(contains(Header, 'AL'));% Columns of the audio loggers
    LogCol = NLogCol(find(cell2mat(DataInfo(NLogCol))==Logger_num)); %#ok<FNDSB>
    if isempty(LogCol) % This is an audiologger and not a neural logger
        LogCol = ALogCol(find(cell2mat(DataInfo(ALogCol))==Logger_num)); %#ok<FNDSB>
        LoggerName{ll} = ['AL' num2str(Logger_num)];
    else
        LoggerName{ll} = ['NL' num2str(Logger_num)];
    end
    BatID{ll} = DataInfo{BatIDCol(find(BatIDCol<LogCol,1,'last'))};
end


% Load data
% Data1 = dir(fullfile(Loggers_dir,['20' Date], sprintf('%s_%s_VocExtractDat*.mat', Date, ExpStartTime)));
% % select the correct files
% Gdf = zeros(length(Data1),1);
% for df=1:length(Data1)
%     if length(strfind(Data1(df).name, '_'))==2
%         Gdf(df)=1;
%     end
% end

% if sum(Gdf)==0
%     warning('No vocalization data extracted by who_calls.m or get_logger_data_voc.m')
% else

%     DataFile = dir(fullfile(Loggers_dir, ['20' Date], sprintf('%s_%s_VocExtractDat*_*.mat', Date, ExpStartTime)));
%     if length(DataFile)~=sum(Gdf)
%         warning('The number of files generated by get_logger_data_voc (Data1) is not the same as the number generated by who_calls (DataFile)')
%         keyboard
%     end
DataFiles = dir(fullfile(Loggers_dir, ['20' Date], sprintf('%s_%s_VocExtractDat*_200.mat', Date, ExpStartTime)));
% Loop through the datafiles
df=0;
while df<length(DataFiles) %1
    df=df+1;
    DataFile = dir(fullfile(Loggers_dir, ['20' Date],sprintf('%s_%s_VocExtractData%d_200.mat', Date, ExpStartTime,df)));
    if isempty(DataFile) % who calls was the earlier format
        DataFile = dir(fullfile(Loggers_dir, ['20' Date], sprintf('%s_%s_VocExtractData_200.mat', Date, ExpStartTime)));
    end
    fprintf(1,'Set %d/%d\nwith file %s\n', df, length(DataFiles), DataFile.name)
    % Save the BatID and Logger names
    Old = load(fullfile(Loggers_dir, ['20' Date], DataFile.name), 'BatID', 'LoggerName');
    if (sum(cell2mat(Old.BatID)==cell2mat(BatID))==length(BatID)) && (sum(strcmp(LoggerName, Old.LoggerName)) == length(BatID))
        fprintf(1,'*** No Error of bat labeling ***\n\n\n')
        % No need to check other files of the same date, the first one
        % is the only one used to initialized these variables
        df=length(DataFiles);
    else
        warning('!!! Error of bat labeling !!!\n')
        %save(fullfile(Loggers_dir, ['20' Date], DataFile.name), 'BatID', 'LoggerName', '-append');

        load(fullfile(Loggers_dir, ['20' Date], DataFile.name),'BioSoundFilenames');
        if ~exist('BioSoundFilenames', 'var')
            fprintf(1,'No vocalization files under biosound\n')
            continue
        end
        %% Loop through calls, correct wavfile and pdf names
        NV = size(BioSoundFilenames,1);
        if isempty(BioSoundFilenames)
            fprintf(1,'No vocalization files under biosound\n')
            continue
        end

        for NVocFile=1:NV
            OldBSName1 = BioSoundFilenames{NVocFile,1};
            OldBSName2=BioSoundFilenames{NVocFile,2};
            OldBSName = OldBSName1;
            if isempty(BioSoundFilenames{NVocFile,1})
                if isempty(BioSoundFilenames{NVocFile,2})
                    fprintf(1,'No vocalization files under biosound\n')
                    continue
                else
                    OldBSName = OldBSName2;
                end
            end
            ALInd = strfind(OldBSName, 'AL');
            ElmtInd = strfind(OldBSName, 'Elmt');
            ALNum = OldBSName((ALInd+2):(ElmtInd-2));

            % ID of the bat
            ALIndex = strcmp(LoggerName, ['AL' ALNum]);
            if sum(ALIndex)~=1
                keyboard
            end
            BatID_local =BatID{ALIndex};
            BatInd = strfind(OldBSName, 'Bat');
            BatID_old = OldBSName((BatInd+3):(ALInd-2));
            if ~strcmp(BatID_old, num2str(BatID_local))
                Error = Error+1;
                %                     if ~isempty(OldBSName1)
                %                         BioSoundFilenames{NVocFile,1}((BatInd+3):(ALInd-2)) = num2str(BatID_local);
                %                         [Status,MSG]=movefile(OldBSName1, BioSoundFilenames{NVocFile,1});
                %                         if ~Status
                %                             warning('Issue with renaming file, error is:')
                %                             MSG %#ok<*NOPRT>
                %                             keyboard
                %                         end
                %                         OldPDFName1 = [OldBSName1(1:end-3) 'pdf'];
                %                         [Status,MSG]=movefile(OldPDFName1, [BioSoundFilenames{NVocFile,1}(1:end-3) 'pdf']);
                %                         if ~Status
                %                             warning('Issue with renaming file, error is:')
                %                             MSG
                %                             keyboard
                %                         end
                %                     end
                %
                %                     BioSoundFilenames{NVocFile,2}((BatInd+3):(ALInd-2)) = num2str(BatID_local);
                %                     [Status,MSG]=movefile(OldBSName2, BioSoundFilenames{NVocFile,2});
                %                     if ~Status
                %                         warning('Issue with renaming file, error is:')
                %                         MSG
                %                         keyboard
                %                     end
                %                     OldPDFName2 = [OldBSName2(1:end-3) 'pdf'];
                %                     [Status,MSG]=movefile(OldPDFName2, [BioSoundFilenames{NVocFile,2}(1:end-3) 'pdf']);
                %                     if ~Status
                %                         warning('Issue with renaming file, error is:')
                %                         MSG
                %                         keyboard
                %                     end
            end
        end
        % save the values!
        %             save(fullfile(Loggers_dir, DataFile.name), 'BioSoundFilenames','-append');
        clear BioSoundFilenames
    end
end
if Error
    fprintf(1,'!!! %d vocalizations are mislabeled with the wrong bat ID !!!\n\n\n', Error)
end
end

function [List2ParamPath] = gather_reconly_datapath(BasePath)
fprintf(1,'*** Gathering paths to audio reconly data ***')
List2ParamPath = cell(10^3,1); % initialize the list to 1000
ExpFolders = dir(fullfile(BasePath,'LMC*'));
NF = 0; % counter for single files
for ee=1:length(ExpFolders)
    fprintf(1, '\n  -> Looking into  %s...\n ', fullfile(ExpFolders(ee).folder,ExpFolders(ee).name))
    DateFolders = dir(fullfile(ExpFolders(ee).folder,ExpFolders(ee).name, 'audio','20*'));
    for dd=1:length(DateFolders)
        fprintf(1, '   %s\n', DateFolders(dd).name);
        AudioParamFiles = dir(fullfile(DateFolders(dd).folder, DateFolders(dd).name,'*RecOnly_param.txt'));
        for ll = 1:length(AudioParamFiles)
            NF = NF +1;
            List2ParamPath{NF} = fullfile(AudioParamFiles(ll).folder, AudioParamFiles(ll).name);
        end
    end
end
List2ParamPath = List2ParamPath(1:NF);
fprintf(1, '\n Files from %d sessions or free recording sessions have been retrieved\n', NF);
end
