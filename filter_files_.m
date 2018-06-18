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
