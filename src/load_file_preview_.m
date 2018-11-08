%--------------------------------------------------------------------------
function mnWav1 = load_file_preview_(fid_bin, P)
    % preview
    if P.nPad_filt > 0
        % [mnWav1, vrWav_mean1, dimm_wav] = load_file_(fid_bin, P.nPad_filt, P);
        dimm_wav = [P.nChans, P.nPad_filt];
        mnWav1 = jrclust.utils.readRecording(fid_bin, P.vcDataType, dimm_wav, ftell(fid_bin), P);
        fseek(fid_bin, -1 * prod(dimm_wav) * bytesPerSample_(P.vcDataType), 'cof');
    else
        mnWav1 = [];
    end
end %func

