%--------------------------------------------------------------------------
function flag = isSingleShank_(P)
    viShank_site = get_(P, 'viShank_site');
    if isempty(viShank_site)
        flag = 1;
    else
        flag = numel(unique(viShank_site)) == 1;
    end
end
