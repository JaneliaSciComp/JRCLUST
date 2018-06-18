%--------------------------------------------------------------------------
% 11/5/17 JJJ: Find .prb file in ./prb/ folder first
% 9/26/17 JJJ: Find .prb file. Look in ./prb/ folder if doesn't exist
function vcFile_prb = find_prb_(vcFile_prb)
    % Find a prb file
    vcFile_prb_ = fullfile(jrcpath_(), 'prb', vcFile_prb);
    if exist_file_(vcFile_prb_)
        vcFile_prb = vcFile_prb_;
    else
        vcFile_prb = search_file_(vcFile_prb, [jrcpath_(), 'prb', filesep()]);
    end
end %func
