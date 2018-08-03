%--------------------------------------------------------------------------
% 12/18/17 JJJ: can handle headerOffset
function nBytes = file_trim_(fid, nBytes, P)
    headerOffset = get_set_(P, 'headerOffset', 0);
    nBytes = nBytes - headerOffset;

    if isempty(P.tlim_load) || ~P.fTranspose_bin
        if headerOffset > 0
            fseek(fid, headerOffset, 'bof');
        end
    else
        bytesPerSample = bytesPerSample_(P.dataType);
        nSamples = floor(nBytes / (bytesPerSample * P.nChans));

        % Apply limit to the range of samples to load
        nlim_load = min(max(round(P.tlim_load * P.sRateHz), 1), nSamples);
        nSamples_load = diff(nlim_load) + 1;
        nBytes = nSamples_load * bytesPerSample * P.nChans;
        fseek_(fid, nlim_load(1), P);
    end
end %func
