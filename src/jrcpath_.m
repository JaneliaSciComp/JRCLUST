%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile_full = jrcpath_(vcFile, fConditional)
    % make it a jrc path
    % Add a full path if the file doesn't exist in the current folder
    %

    if nargin<1, vcFile = ''; end
    if nargin<2, fConditional = 0; end

    jrcpath = fileparts(mfilename('fullpath'));
    if fConditional
        if exist(vcFile, 'file') == 2
            vcFile_full = vcFile;
            return;
        end
    end
    vcFile_full = fullfile(jrcpath, vcFile);
    % if exist(vcFile_full, 'file') ~= 2
    %     vcFile_full = [];
    % end
end %func
