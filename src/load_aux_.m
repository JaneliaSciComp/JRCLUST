function [auxTraces, auxTime] = load_aux_(hCfg)
    [auxTraces, auxTime] = deal([]);
    if numel(hCfg.rawRecordings) > 1
        jrclust.utils.qMsgBox('Multi-file mode is currently not supported');
        return;
    end

    % try to guess auxFile
    if isempty(hCfg.auxFile)
        [~, ~, ext] = fileparts(hCfg.rawRecordings{1});

        if strcmpi(ext, '.ns5')
            try
                hCfg.auxFile = jrclust.utils.subsExt(hCfg.rawRecordings{1}, '.ns2');
            catch ME
                fprintf(2, '%s\n', ME.message);
                return;
            end
        elseif ismember(lower(ext), {'.bin', '.dat'})
            try
                hCfg.auxFile = hCfg.rawRecordings{1};
            catch ME
                fprintf(2, '%s\n', ME.message);
                return;
            end
        else
            fprintf(2, 'Unable to determine the aux file. Please set "auxFile" in %s.\n', hCfg.configFile);
            return;
        end
    end

    [~, ~, auxExt] = fileparts(hCfg.auxFile);
    switch lower(auxExt)
%         case '.ns2'
%             auxChan = hCfg.getOr('auxChan', 1);
%             [mnWav_aux, hFile_aux, auxData] = load_nsx_(hCfg.auxFile);
%             scale_aux = hFile_aux.Entity(auxChan).Scale * hCfg.auxScale;
%             vrWav_aux = single(mnWav_aux(auxChan,:)') * scale_aux;
%             auxRate = auxData.sRateHz;
        case '.mat'
            auxData = load(hCfg.auxFile);
            auxDataFields = fieldnames(auxData);
            auxTraces = auxData.(auxDataFields{1});
            auxRate = hCfg.getOr('auxRate', hCfg.sampleRate);
        case {'.dat', '.bin'}
            if isempty(hCfg.auxChan)
                return;
            end

            hRec = jrclust.models.recording.Recording(hCfg.auxFile, hCfg);
            auxTraces = single(hRec.readROI(hCfg.auxChan, 1:hRec.nSamples))*hCfg.bitScaling*hCfg.auxScale;
            auxRate = hCfg.getOr('auxRate', hCfg.sampleRate);
        otherwise
            fprintf(2, 'hCfg.auxFile: unsupported file format: %s\n', auxExt);
        return;
    end % switch

    if nargout >= 2
        auxTime = single(1:numel(auxTraces))' / auxRate;
    end
end %func
