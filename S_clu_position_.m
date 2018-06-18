%--------------------------------------------------------------------------
function S_clu = S_clu_position_(S_clu, viClu_update)
    % determine cluster position from spike position
    % 6/27/17 JJJ: multiple features supported (single dimension such as energy and Vpp)
    global trFet_spk
    if nargin<2, viClu_update = []; end
    P = get0_('P'); %P = S_clu.P;
    if ~isfield(S_clu, 'vrPosX_clu'), S_clu.vrPosX_clu = []; end
    if ~isfield(S_clu, 'vrPosY_clu'), S_clu.vrPosY_clu = []; end

    if isempty(S_clu.vrPosX_clu) || ~isempty(S_clu.vrPosY_clu)
        viClu_update = [];
    end
    if isempty(viClu_update)
        [vrPosX_clu, vrPosY_clu] = deal(zeros(S_clu.nClu, 1));
        viClu1 = 1:S_clu.nClu;
    else % selective update
        vrPosX_clu = S_clu.vrPosX_clu;
        vrPosY_clu = S_clu.vrPosY_clu;
        viClu1 = viClu_update(:)';
    end
    viSites_fet = 1:(1+P.maxSite*2-P.nSites_ref);
    for iClu = viClu1
        %     viSpk_clu1 = S_clu.cviSpk_clu{iClu};
        [viSpk_clu1, viSites_clu1] = S_clu_subsample_spk_(S_clu, iClu);
        if isempty(viSpk_clu1), continue; end

        viSites_clu1 = viSites_clu1(1:end-P.nSites_ref);
        mrVp1 = squeeze_(trFet_spk(viSites_fet,1,viSpk_clu1));
        mrSiteXY1 = single(P.mrSiteXY(viSites_clu1,:)); %electrode

        vrPosX_clu(iClu) = median(centroid_mr_(mrVp1, mrSiteXY1(:,1), 2));
        vrPosY_clu(iClu) = median(centroid_mr_(mrVp1, mrSiteXY1(:,2), 2));
    end
    S_clu.vrPosX_clu = vrPosX_clu;
    S_clu.vrPosY_clu = vrPosY_clu;
end %func
