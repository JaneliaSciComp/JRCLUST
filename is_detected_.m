%--------------------------------------------------------------------------
function flag = is_detected_(P)
    % return true if already detected. .spkwav file must exist
    vcFile = strrep(P.paramFile, '.prm', '_waveforms.bin');
    flag = fileExists(vcFile);
    if flag, flag = getBytes_(vcFile) > 0; end
end %func
