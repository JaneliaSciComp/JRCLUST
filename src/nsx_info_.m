%--------------------------------------------------------------------------
% 9/22/17 JJJ: Created for SPARC
function [P, nSamples, hFile] = nsx_info_(vcFile_nsx)
    addpath('./neuroshare/');
    [ns_RESULT, hFile] = ns_OpenFile(vcFile_nsx);
    vlAnalog_chan= strcmpi({hFile.Entity.EntityType}, 'Analog');
    nSamples = hFile.TimeSpan / hFile.FileInfo.Period;
    % viElecID = double([hFile.Entity.ElectrodeID]);
    P = struct('vcDataType', 'int16', 'nChans', sum(vlAnalog_chan), ...
    'uV_per_bit', hFile.Entity(1).Scale, 'sRateHz', 30000 / hFile.FileInfo.Period);
end %func
