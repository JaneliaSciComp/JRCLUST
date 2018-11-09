%--------------------------------------------------------------------------
% 8/4/17 JJJ: Accepts a text file containing a list of bin files
% 7/31/17 JJJ: handles wild card combined with cells
function [csFiles_valid, viValid] = filter_files_(csFiles, fSortMode)
    % Files must exist and non-zero byte
    if isTextFile_(csFiles)
        csFiles = load_batch_(csFiles);
    else
        if nargin<2
            P = get0_('P');
            fSortMode = get_set_(P, 'sort_file_merge', 1);
        end
        csFiles = list_files_(csFiles, fSortMode);
    end

    % filter based on file presence and bytes
    vlValid = false(size(csFiles));
    for iFile=1:numel(csFiles)
        S_dir1 = dir(csFiles{iFile});
        if isempty(S_dir1), continue; end
        if S_dir1.bytes == 0, continue; end
        vlValid(iFile) = 1;
    end
    viValid = find(vlValid);
    csFiles_valid = csFiles(viValid);
    if ~all(vlValid)
        fprintf('Files not found:\n');
        cellfun(@(vc)fprintf('\t%s\n', vc), csFiles(~vlValid));
    end
end %func

%% local functions
%--------------------------------------------------------------------------
% 7/31/17 JJJ: Expand cellstring and wild card to list of files
function csFiles = list_files_(csFiles, fSortMode)
    if nargin < 2
        P = get0_('P');
        fSortMode = get_set_(P, 'sort_file_merge', 1);
    end
    if ischar(csFiles)
        if any(csFiles == '*')
            csFiles = dir_file_(csFiles, fSortMode); %sort by dates
        else
            csFiles = {csFiles}; %make it a cell
        end
    elseif iscell(csFiles)
        csFiles = cellfun(@(vc)dir_file_(vc, fSortMode), csFiles, 'UniformOutput', 0);
        csFiles = [csFiles{:}];
    end
end %func
