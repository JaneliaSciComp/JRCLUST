%--------------------------------------------------------------------------
function varargout = trWav2fet_cov_(trWav2, P)
    % tnWav1: nT x nSites_spk x nSpk
    % subtract ref
    % nSites_spk = 1 + 2 * P.maxSite - P.nSites_ref; % size(tnWav_spk, 2);
    vnDelay_fet = get_set_(P, 'vnDelay_fet', [0,3]);

    [nT, nSpk, nSites_spk] = size(trWav2);
    [cvi2, cvi1] = shift_range_(nT, [], vnDelay_fet);
    cmrFet = cell(numel(vnDelay_fet), 1);
    for iDelay = 1:numel(vnDelay_fet)
        mr1_ = meanSubt_(trWav2(cvi1{iDelay},:,1));
        mr1_ = bsxfun(@rdivide, mr1_, sqrt(mean(mr1_.^2))); %zscore fast
        tr1_ = repmat(mr1_, [1,1,nSites_spk]);
        tr2_ = meanSubt_(trWav2(cvi2{iDelay},:,:));
        cmrFet{iDelay} = permute(mean(tr1_ .* tr2_, 1), [3,2,1]);
    end %for

    if nargout==1
        varargout{1} = cell2mat_(cmrFet);
    else
        varargout = cmrFet;
    end
end %func
