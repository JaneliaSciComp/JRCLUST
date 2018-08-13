%--------------------------------------------------------------------------
% 9/27/17 Validate parameters
function flag = validate_param_(P)
    % validate P

    NDIM_SORT_MAX = get_set_(P, 'nC_max', 45);

    csError = {};

    nSites_spk = P.maxSite * 2 + 1 - P.nSites_ref;
    nFet_sort = nSites_spk * P.nPcPerChan;

    if nSites_spk <= 0
        csError{end+1} = sprintf('Negative # Sites/spk. Use this formula to adjust maxSite and nSites_ref (nSites_spk = 1+maxSite*2-nSites_ref)');
    end
    if nFet_sort > NDIM_SORT_MAX
        csError{end+1} = sprintf('# dimensions (%d) exceeds the maximum limit for CUDA code (%d), decrease maxSite', nFet_sort, NDIM_SORT_MAX);
    end

    % Validate format


    % Validate display


    % Validate preprocess


    % Validate feature


    % Validate cluster


    % Validate post-cluster


    if isempty(csError)
        flag = true;
    else
        cellfun(@(vc_)fprintf(2, '%s\n', vc_), csError);
        flag = false;
    end
end %func
