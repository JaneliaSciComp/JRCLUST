%--------------------------------------------------------------------------
function [viSite, viSite2, viSite3] = S_clu_peak2_(S_clu)
    mrMin_clu = squeeze_(min(S_clu.tmrWav_spk_clu) - S_clu.tmrWav_spk_clu(1,:,:));
    % mrMin_clu = squeeze_(min(S_clu.tmrWav_spk_clu) - max(S_clu.tmrWav_spk_clu));
    % mrMin_clu = squeeze_(min(S_clu.tmrWav_spk_clu));
    % mrMin_clu = squeeze_(min(S_clu.tmrWav_raw_clu));

    [~, viSite] = min(mrMin_clu);
    if nargout>=2
        mrMin_clu(sub2ind(size(mrMin_clu), viSite, 1:numel(viSite))) = 0;
        [~, viSite2] = min(mrMin_clu);
    end
    if nargout>=3
        mrMin_clu(sub2ind(size(mrMin_clu), viSite2, 1:numel(viSite2))) = 0;
        [~, viSite3] = min(mrMin_clu);
    end
end % function
