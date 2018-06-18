%--------------------------------------------------------------------------
% 9/26/17 JJJ: Merged to Master. This function is not used, instead JRCLUST deals with the offset
% 9/22/17 JJJ: Created for SPARC.
function [P, nSamples, vcFile_bin] = nsx2bin_(vcFile_nsx, fInvert)
    if nargin<2, fInvert = 0; end

    nBuffer = 1e8; % in bytes
    addpath('./neuroshare/');
    vcFile_bin = subsFileExt_(vcFile_nsx, '.bin');
    [ns_RESULT, hFile] = ns_OpenFile(vcFile_nsx);
    % [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    vlAnalog_chan= strcmpi({hFile.Entity.EntityType}, 'Analog');
    nSamples = hFile.TimeSpan / hFile.FileInfo.Period;
    nChans = sum(vlAnalog_chan);
    % viElecID = double([hFile.Entity.ElectrodeID]);
    P = struct('vcDataType', 'int16', 'nChans', nChans, ...
    'uV_per_bit', hFile.Entity(1).Scale, 'sRateHz', 30000 / hFile.FileInfo.Period);
    fprintf('Loading %s...', vcFile_nsx); t_load = tic;

    fid = hFile.FileInfo.FileID;
    fidw = fopen(vcFile_bin, 'w');

    if ismember(hFile.FileInfo.FileTypeID, {'NEURALCD', 'NEUCDFLT'})
        fseek(fid, hFile.FileInfo.BytesHeaders + 9, -1);
    else
        fseek(fid, hFile.FileInfo.BytesHeaders, -1);
    end
    nLoad = ceil(nSamples*nChans / nBuffer);
    for iLoad = 1:nLoad
        if iLoad==nLoad
            nBuffer_ = nSamples*nChans - (nLoad-1) * nBuffer;
        else
            nBuffer_ = nBuffer;
        end
        vnBuffer = fread(fid, nBuffer_, '*int16');
        if fInvert, vnBuffer = -vnBuffer; end
        fwrite(fidw, vnBuffer, 'int16');
    end
    fclose(fid);
    fclose(fidw);
    fprintf('took %0.1fs\n', toc(t_load));
end %func
