%--------------------------------------------------------------------------
% 11/5/17 JJJ: Find .prb file in ./prb/ folder first
% 9/26/17 JJJ: Find .prb file. Look in ./prb/ folder if doesn't exist
function vcFile_prb = find_prb_(vcFile_prb)
    if nargin < 1 || isempty(vcFile_prb)
        error('Probe file not given');
    end

    % check the probe directory
    probedir = fullfile(fileparts(jrcpath_()), 'prb');
    vcFile_prb = search_file_(vcFile_prb, probedir);
end %func
