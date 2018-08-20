%--------------------------------------------------------------------------
function fseek_(fid_bin, iSample_bin, P)
    % provide # of samples to skip for transpose multi-channel type
    % Index starts with 1 unlike fseek() function

    if ~P.fTranspose_bin, return; end
    headerOffset = getOr(P, 'headerOffset', 0);
    iOffset = (iSample_bin-1) * P.nChans * bytesPerSample_(P.dataType) + headerOffset;
    if iOffset<0, iOffset = 0; end
    fseek(fid_bin, iOffset, 'bof');
end %func
