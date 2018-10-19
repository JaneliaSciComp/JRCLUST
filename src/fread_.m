%--------------------------------------------------------------------------
function mnWav1 = fread_(fid_bin, dimm_wav, vcDataType)
    % Get around fread bug (matlab) where built-in fread resize doesn't work
    try
        if isempty(dimm_wav)
            mnWav1 = fread(fid_bin, inf, ['*', vcDataType]);
        else
            if numel(dimm_wav)==1, dimm_wav = [dimm_wav, 1]; end
            mnWav1 = fread(fid_bin, prod(dimm_wav), ['*', vcDataType]);
            if numel(mnWav1) == prod(dimm_wav)
                mnWav1 = reshape(mnWav1, dimm_wav);
            else
                dimm2 = floor(numel(mnWav1) / dimm_wav(1));
                if dimm2 >= 1
                    mnWav1 = reshape(mnWav1, dimm_wav(1), dimm2);
                else
                    mnWav1 = [];
                end
            end
        end
    catch
        disperr_();
    end
end %func
