%--------------------------------------------------------------------------
function S_clu = S_clu_cleanup_(S_clu, P)
    % 17/7/3: Cluster cleanup routine, Mahal distance based outlier removal

    viSite2_spk = get0_('viSite2_spk');
    thresh_mad_clu = get_set_(P, 'thresh_mad_clu', 7.5);
    if thresh_mad_clu == 0, return; end % aborted

    fprintf('Cleaning up clusters\n\t'); t1=tic;
    warning off;
    for iClu = 1:S_clu.nClu
        mrFet1_clu1 = S_clu_getFet_(S_clu, iClu, viSite2_spk)';
        viSpk_clu1 = S_clu.cviSpk_clu{iClu};
        try
            vrDist_clu1 = madscore_(log(mahal(mrFet1_clu1, mrFet1_clu1)));
        catch
            continue;
        end
        vlExcl_clu1 = (vrDist_clu1 > thresh_mad_clu);
        if any(vlExcl_clu1)
            S_clu.cviSpk_clu{iClu} = viSpk_clu1(~vlExcl_clu1);
            S_clu.viClu(viSpk_clu1(vlExcl_clu1)) = 0; %classify as noise cluster
            S_clu.vnSpk_clu(iClu) = numel(S_clu.cviSpk_clu{iClu});
        end
        fprintf('.');
    end %for
    fprintf('\n\ttook %0.1fs.\n', toc(t1));
end %func
