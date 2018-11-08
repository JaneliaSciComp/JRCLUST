%--------------------------------------------------------------------------
function mnWav = load_bin_multi_(fid, cvi_lim_bin, P)
    mnWav = cell(size(cvi_lim_bin));
    fpos = 0;
    for i=1:numel(cvi_lim_bin)
        lim1 = cvi_lim_bin{i};
        if i>1, fseek_(fid, lim1(1), P); end
        mnWav{i} = jrclust.utils.readBin(fid, P.vcDataType, [P.nChans, diff(lim1)+1]);
        if i==1, fpos = ftell(fid); end
    end %for
    mnWav = cell2mat_(mnWav);
    fseek(fid, fpos, 'bof'); % restore the file position
end %func
