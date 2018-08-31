%--------------------------------------------------------------------------
% 8/15/18 ACL: Check if is an absolute path before anything else
% 11/5/17 JJJ: Find .prb file in ./prb/ folder first
% 9/26/17 JJJ: Find .prb file. Look in ./prb/ folder if doesn't exist
function probeFilename = find_prb_(probeFilename)
    % Find a prb file
    if ~isAbsPath(probeFilename)
        vcFile_prb_ = fullfile(jrcpath_(), 'prb', probeFilename);

        if fileExists(vcFile_prb_)
            probeFilename = vcFile_prb_;
        else
            probeFilename = search_file_(probeFilename, [jrcpath_(), 'prb', filesep()]);
        end
    end
end % function
