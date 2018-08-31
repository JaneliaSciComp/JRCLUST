%--------------------------------------------------------------------------
function S_clu = S_clu_quality_(S_clu, P, viClu_update)
    % 7/5/17 JJJ: Added isolation distance, L-ratio,
    % TODO: update when deleting, splitting, or merging
    if nargin<3, viClu_update = []; end
    t1 = tic;
    fprintf('Calculating cluster quality...\n');
    % [vrVmin_clu, vrVmax_clu, vrVpp_clu] = clu_amp_(S_clu.tmrWav_clu, S_clu.clusterSites);
    % [vrVmin_uV_clu, vrVmax_uV_clu, vrVpp_uV_clu] = clu_amp_(S_clu.tmrWav_raw_clu, S_clu.clusterSites);
    mrVmin_clu = squeeze_(min(S_clu.tmrWav_clu));
    mrVmax_clu = squeeze_(max(S_clu.tmrWav_clu));
    mrVmin_uv_clu = squeeze_(min(S_clu.tmrWav_raw_clu));
    mrVmax_uv_clu = squeeze_(max(S_clu.tmrWav_raw_clu));
    % tmrWav_clu = S_clu.tmrWav_clu;
    % mrVmin_clu = shiftdim(min(tmrWav_clu,[],1));
    % mrVmax_clu = shiftdim(max(tmrWav_clu,[],1));
    vrVpp_clu = mr2vr_sub2ind_(mrVmax_clu-mrVmin_clu, S_clu.clusterSites, 1:S_clu.nClusters);
    vrVmin_clu = mr2vr_sub2ind_(mrVmin_clu, S_clu.clusterSites, 1:S_clu.nClusters);
    vrVpp_uv_clu = mr2vr_sub2ind_(mrVmax_uv_clu-mrVmin_uv_clu, S_clu.clusterSites, 1:S_clu.nClusters);
    vrVmin_uv_clu = mr2vr_sub2ind_(mrVmin_uv_clu, S_clu.clusterSites, 1:S_clu.nClusters);
    % [vrVpp_clu, ~] = max(mrVmax_clu - mrVmin_clu,[],1);
    % [vrVmin_clu, viSite_min_clu] = min(mrVmin_clu,[],1);
    % if ~isfield(S_clu, 'clusterSites')
    %     clusterSites = viSite_min_clu;
    % else
    %     clusterSites = S_clu.clusterSites;
    % end
    % vrVpp_clu=vrVpp_clu(:);  clusterSites=clusterSites(:);
    % [vrVpp_clu, clusterSites, vrVmin_clu] = multifun_(@(x)x(:), vrVpp_clu, clusterSites, vrVmin_clu);

    try
        S0 = get(0, 'UserData');
        if isempty(S0), S0 = load0_(subsFileExt_(P.paramFile, '_jrc.mat')); end
        vrVrms_site = bit2uV_(single(S0.vrThresh_site(:)) / S0.P.qqFactor, P);
        %     vrSnr_clu = vrVpp_clu ./ vrVrms_site(clusterSites);
        vrSnr_clu = abs(vrVmin_clu) ./ vrVrms_site(S_clu.clusterSites);
        vnSite_clu = sum(bsxfun(@lt, mrVmin_clu, -vrVrms_site * S0.P.qqFactor),1)';
    catch
        [vrVrms_site, vrSnr_clu, vnSite_clu] = deal([]);
        disp('no Sevt in memory.');
    end
    [vrIsoDist_clu, vrLRatio_clu, vrIsiRatio_clu] = S_clu_quality2_(S_clu, P, viClu_update);

    S_clu = struct_add_(S_clu, vrVpp_clu, vrSnr_clu, vrVrms_site, vnSite_clu, ...
    vrIsoDist_clu, vrLRatio_clu, vrIsiRatio_clu, vrVpp_uv_clu, vrVmin_uv_clu);
    fprintf('\ttook %0.1fs\n', toc(t1));
end % function
