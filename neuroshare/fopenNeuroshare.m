%--------------------------------------------------------------------------
% 9/22/17 JJJ: Created for SPARC
function [fid, nBytes, headerOffset] = fopenNeuroshare(vcFile_nsx)
    [ns_RESULT, hFile] = ns_OpenFile(vcFile_nsx);
    % [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    vlAnalog_chan= strcmpi({hFile.Entity.EntityType}, 'Analog');
    nSamples = hFile.TimeSpan / hFile.FileInfo.Period;
    nChans = sum(vlAnalog_chan);
    nBytes = nSamples * nChans * 2;
    fid = hFile.FileInfo.FileID;

    if ismember(hFile.FileInfo.FileTypeID, {'NEURALCD', 'NEUCDFLT'})
        headerOffset = hFile.FileInfo.BytesHeaders + 9;
    else
        headerOffset = hFile.FileInfo.BytesHeaders;
    end
    fseek(fid, headerOffset, -1);
end % function
