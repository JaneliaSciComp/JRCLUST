%--------------------------------------------------------------------------
function fseek_(fid_bin, iSample_bin, P)
    % provide # of samples to skip for transpose multi-channel type
    % Index starts with 1 unlike fseek() function

    if P.fTranspose_bin
        fseek(fid_bin, max(0, (iSample_bin-1) * P.nChans * bytesPerSample_(P.vcDataType) + get_set_(P, 'header_offset', 0)), 'bof');
    end
end %func
