function P = doCreateConfigFile(outputDir, binFiles, probeFile, templateFile, fAsk)
    %DOCREATECONFIGFILE
    templateFile = jrclust.utils.absPath(templateFile);
    if ~isempty(templateFile)
        P0 = jrclust.utils.mToStruct(templateFile); % a struct on success or [] on failure
    else
        P0 = [];
    end

    P.rawRecordings = binFiles;

    % load meta file
    metaFiles = cellfun(@(bf) jrclust.utils.subsExt(bf, '.meta'), binFiles, 'UniformOutput', 0);
    metaFiles = metaFiles(cellfun(@(fn) exist(fn, 'file') == 2, metaFiles));
    if ~isempty(metaFiles)
        metaFile = metaFiles{1}; % just take the first one
        SMeta = jrclust.utils.loadMetadata(metaFile);
    else
        SMeta = [];
    end

    % assign config file name
    [~, probeName, ~] = fileparts(probeFile);
    [~, sessionName, ~] = fileparts(binFiles{1});
    sessionName = [sessionName '_', probeName];
    P.configFile = fullfile(outputDir, [sessionName '.prm']);

    % set probe file
    P.probeFile = probeFile;

    if exist(P.configFile, 'file') && fAsk
        dlgAns = questdlg('File already exists. Overwrite prm file?', 'Warning', 'Yes', 'No', 'No');
        if ~strcmp(dlgAns, 'Yes')
            P = [];
            return;
        end
    elseif ~exist(P.configFile, 'file') % create it
        fid = fopen(P.configFile, 'wt');
        if fid == -1
            error('Could not create config file %s', P.configFile);
        end

        fclose(fid);
    end

    P = jrclust.utils.mergeStructs(P0, P);
    P = jrclust.utils.mergeStructs(P, SMeta);
end
