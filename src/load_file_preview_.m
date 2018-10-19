%--------------------------------------------------------------------------
function mnWav1 = load_file_preview_(fid_bin, P)
    % preview
    if P.nPad_filt > 0
        [mnWav1, vrWav_mean1, dimm_wav] = load_file_(fid_bin, P.nPad_filt, P);
        frewind_(fid_bin, dimm_wav, P.vcDataType);
    else
        mnWav1 = [];
    end
end %func
