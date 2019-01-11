function [auxSamples, auxTimes] = loadAuxChannel(hCfg)
    %LOADAUXCHANNEL Load the aux channel 
    [auxSamples, auxTimes] = deal([]);
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
            catch
                return;
            end
        elseif ismember(lower(ext), {'.bin', '.dat'})
            try
                hCfg.auxFile = hCfg.rawRecordings{1};
            catch
                return;
            end
        else
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
            auxSamples = auxData.(auxDataFields{1});
            auxRate = hCfg.getOr('auxRate', hCfg.sampleRate);
        case {'.dat', '.bin'}
            if isempty(hCfg.auxChan)
                return;
            end

            hRec = jrclust.models.recording.Recording(hCfg.auxFile, hCfg);
            auxSamples = single(hRec.readROI(hCfg.auxChan, 1:hRec.nSamples))*hCfg.bitScaling*hCfg.auxScale;
            auxRate = hCfg.getOr('auxRate', hCfg.sampleRate);
        otherwise
            jrclust.utils.qMsgBox(sprintf('hCfg.auxFile: unsupported file format: %s\n', auxExt));
        return;
    end % switch

    if nargout >= 2
        auxTimes = single(1:numel(auxSamples))'/auxRate;
    end
end
