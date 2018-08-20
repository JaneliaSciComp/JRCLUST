%--------------------------------------------------------------------------
% 9/19/17 JJJ: Created for SPARC
function [mnWav, hFile, P] = load_nsx_(vcFile_nsx)
    addpath('./neuroshare/');
    [ns_RESULT, hFile] = ns_OpenFile(vcFile_nsx);
    % [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    vlAnalog_chan= strcmpi({hFile.Entity.EntityType}, 'Analog');
    nSamples = hFile.TimeSpan / hFile.FileInfo.Period;
    nChans = sum(vlAnalog_chan);
    % viElecID = double([hFile.Entity.ElectrodeID]);
    P = struct('dataType', 'int16', 'nChans', nChans, ...
    'uV_per_bit', hFile.Entity(1).Scale, 'sampleRateHz', 30000 / hFile.FileInfo.Period);

    fprintf('Loading %s...', vcFile_nsx); t_load = tic;
    fid = hFile.FileInfo.FileID;
    if ismember(hFile.FileInfo.FileTypeID, {'NEURALCD', 'NEUCDFLT'})
        fseek(fid, hFile.FileInfo.BytesHeaders + 9, -1);
    else
        fseek(fid, hFile.FileInfo.BytesHeaders, -1);
    end
    mnWav = fread(fid, [nChans, nSamples], '*int16');
    fclose(fid);
    fprintf('took %0.1fs\n', toc(t_load));

    if nargout==0, assignWorkspace_(mnWav, hFile); end
end % func
