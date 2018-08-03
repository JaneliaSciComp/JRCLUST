%--------------------------------------------------------------------------
function frewind_(fid_bin, dimm_wav, dataType)
    % move the file pointer back by the dimm_wav, vcDatatype
    fseek(fid_bin, -1 * prod(dimm_wav) * bytesPerSample_(dataType), 'cof');
end %func
