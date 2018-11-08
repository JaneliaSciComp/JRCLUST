%--------------------------------------------------------------------------
function bin_file_copy_(vcFile_r, vcFile_w, P)
    % bin_file_copy_(vcFile_r, fid_w, P): provide file handle to write
    % bin_file_copy_(vcFile_r, vcFile_w, P): provide file name to write

    if ischar(vcFile_w)
        fid_w = fopen(vcFile_w, 'w');
    else
        fid_w = vcFile_w;
    end
    try
        vn = jrclust.utils.readBin(vcFile_r, P.vcDataType, [], P.header_offset);
        if P.nSkip_lfp > 1
            nSkip_lfp = round(P.nSkip_lfp);
            vn = reshape_(vn, P.nChans)'; % return nTime x nChans
            dimm_ = size(vn); % [nTime x nChans]
            nSamples_lfp = floor(dimm_(1) / nSkip_lfp);
            if dimm_(1) ~= (nSamples_lfp * nSkip_lfp)
                vn = vn(1:nSamples_lfp * nSkip_lfp, :); % trim extra time
            end
            vn = mean(reshape(vn, nSkip_lfp, []), 1);
            vn = reshape(cast(vn, P.vcDataType), nSamples_lfp, [])';
        end
        write_bin_(fid_w, vn); %write all straight
    catch E
        disperr_('bin_file_copy_', E);
    end
    if ischar(vcFile_w)
        fclose(fid_w);
    end
end %func
