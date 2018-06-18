%--------------------------------------------------------------------------
% 10/10/17 JJJ: created and tested
function tnWav_spk1 = fread_spk_(fid_, dimm_, viSpk1);
    % ~35x slower than RAM indexing

    tnWav_spk1 = zeros([dimm_(1), dimm_(2), numel(viSpk1)], 'int16');
    fseek(fid_, 0, 'bof');
    vnOffset = (diff(viSpk1) - 1) * dimm_(1) * dimm_(2) * 2;
    for iSpk = 1:numel(viSpk1)
        if iSpk>1, fseek(fid_, vnOffset(iSpk-1), 'cof'); end
        tnWav_spk1(:,:,iSpk) = fread(fid_, dimm_(1:2), '*int16');
    end
end %func
