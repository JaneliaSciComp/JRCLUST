%--------------------------------------------------------------------------
% 12/18/17 JJJ: can handle header_offset
function nBytes_load = file_trim_(fid, nBytes_load, P)
    header_offset = get_set_(P, 'header_offset', 0);
    nBytes_load = nBytes_load - header_offset;
    if isempty(P.tlim_load) || ~P.fTranspose_bin
        if header_offset > 0, fseek(fid, header_offset, 'bof'); end
        return;
    end

    bytesPerSample = bytesPerSample_(P.vcDataType);
    nSamples = floor(nBytes_load / bytesPerSample / P.nChans);

    % Apply limit to the range of samples to load
    nlim_load = min(max(round(P.tlim_load * P.sRateHz), 1), nSamples);
    nSamples_load = diff(nlim_load) + 1;
    nBytes_load = nSamples_load * bytesPerSample * P.nChans;
    % if nlim_load(1)>1,
    fseek_(fid, nlim_load(1), P);
    % end
end %func
