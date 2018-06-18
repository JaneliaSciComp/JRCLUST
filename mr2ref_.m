%--------------------------------------------------------------------------
% 8/17/17 JJJ: created
function vrRef = mr2ref_(mnWav_filt, vcCommonRef, viSite_bad)
    if nargin<2, vcCommonRef = 'mean'; end
    if nargin<3, viSite_bad = []; end
    if ~isempty(viSite_bad)
        mnWav_filt(:,viSite_bad) = [];
    end
    if strcmpi(vcCommonRef, 'median')
        vrRef = median(mnWav_filt, 2);
    else
        vrRef = mean(mnWav_filt, 2);
    end
end %func
