function manualType_calls(Loggers_dir, Date, ExpStartTime)
% This function is run on all data to check if call are Trills or none
% trills
% Load data
DataFiles = dir(fullfile(Loggers_dir, sprintf('%s_%s_VocExtractDat*_*.mat', Date, ExpStartTime)));

% Loop through the datafiles
for df=1:length(DataFiles) %1
    fprintf(1, 'Set %d/%d\n', df, length(DataFiles))
    DataFile = DataFiles(df);
    load(fullfile(DataFile.folder, DataFile.name), 'ManualCallType');
    if exist('ManualCallType', 'var') && (isempty(ManualCallType) || ~isempty(ManualCallType(end)))
        clear ManualCallType
        continue
        
    elseif exist('ManualCallType', 'var') && isempty(ManualCallType(end))
        GoManualCallType = input(sprintf('It looks like we should start from here because the last vocalization is labelled as NaN in ManualCallTYpe\n There is however %d Nan\n Resume ManualCallType:1 skip and check the next file:0\n', sum(isnan(ManualCallType))));
        if GoManualCallType
            Rangevv = find(isnan(ManualCallType));
            load(fullfile(DataFile.folder, DataFile.name), 'BioSoundCalls')
            NVoc = size(BioSoundCalls,1);
        else
            clear ManualCallType
            continue
        end
    else
        load(fullfile(DataFile.folder, DataFile.name), 'BioSoundCalls')
        if ~exist('BioSoundCalls', 'var') || isempty(BioSoundCalls)
            fprintf(1, 'No Calls for that set\n')
            ManualCallType = [];
            save(fullfile(DataFile.folder, DataFile.name), 'ManualCallType', '-append')
            clear ManualCallType BioSoundCalls
            continue
        end
        NVoc = size(BioSoundCalls,1);
        ManualCallType = cell(NVoc,1);
        Rangevv = 1:NVoc;
    end
    for jj=1:length(Rangevv)
        vv=Rangevv(jj);
        fprintf('set %d/%d Voc %d/%d:  ',df,length(DataFiles), vv,NVoc)
        if isempty(BioSoundCalls{vv,1})
            warning('Uncorrect sound extract that has no Biosound values')
            continue
        else
            APM=audioplayer(BioSoundCalls{vv,1}.sound./(max(abs(BioSoundCalls{vv,1}.sound))),BioSoundCalls{vv,1}.samprate);
            APP=audioplayer(BioSoundCalls{vv,2}.sound./(max(abs(BioSoundCalls{vv,2}.sound))),BioSoundCalls{vv,2}.samprate);
            
%             pause(1)
%             AP2=audioplayer(BioSoundCalls{vv,2}.sound./(max(abs(BioSoundCalls{vv,2}.sound))),BioSoundCalls{vv,2}.samprate);
%             play(AP2)
            INPUT=[];
            while isempty(INPUT)
                play(APM)
                pause(1)
                play(APP)
                INPUT = input('Trill (1), Bark(2), pitchy call(3), low buzz(4) or Unknown (0)');
                if isempty(INPUT) || ((INPUT~=0) &&(INPUT~=1) && (INPUT~=2) && (INPUT~=3)&& (INPUT~=4))
                    INPUT=[];
                end
            end
            if INPUT==1
                ManualCallType{vv} = 'Tr';
            elseif INPUT==2
                ManualCallType{vv} = 'Ba';
            elseif INPUT==3
                ManualCallType{vv} = 'Pi';
            elseif INPUT==4
                ManualCallType{vv} = 'Bu';
            elseif INPUT==0
                ManualCallType{vv} = 'Un';
            end
        end
%             keyboard
    end
    save(fullfile(DataFile.folder, DataFile.name), 'ManualCallType', '-append')
    
    clear BioSoundCalls ManualCallType
end
end