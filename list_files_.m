%--------------------------------------------------------------------------
% 7/31/17 JJJ: Expand cellstring and wild card to list of files
function csFiles = list_files_(csFiles, fSortMode)
    if nargin<2
        P = get0_('P');
        fSortMode = getOr(P, 'sort_file_merge', 1);
    end
    if ischar(csFiles)
        if any(csFiles=='*')
            csFiles = dir_file_(csFiles, fSortMode); %sort by dates
        else
            csFiles = {csFiles}; %make it a cell
        end
    elseif iscell(csFiles)
        csFiles = cellfun(@(vc)dir_file_(vc, fSortMode), csFiles, 'UniformOutput', 0);
        csFiles = [csFiles{:}];
    end
end %func
