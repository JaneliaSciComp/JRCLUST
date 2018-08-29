%--------------------------------------------------------------------------
function vi = cell2vi_(cvi)
    % convert cell index to array of index
    vn_site = cellfun(@(x)numel(x), cvi); %create uniform output
    vi = cell(numel(cvi), 1);
    for iSite=1:numel(cvi)
        vi{iSite} = iSite * ones(vn_site(iSite), 1);
    end
    vi = cell2mat_(vi);
end % function
