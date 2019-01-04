%--------------------------------------------------------------------------
function S_clu = S_clu_quality_(S_clu, P, viClu_update)
    % 7/5/17 JJJ: Added isolation distance, L-ratio,
    % TODO: update when deleting, splitting, or merging
    if nargin<3, viClu_update = []; end
    t1 = tic;
    fprintf('Calculating cluster quality...\n');
    % [vrVmin_clu, vrVmax_clu, vrVpp_clu] = clu_amp_(S_clu.tmrWav_clu, S_clu.viSite_clu);
    % [vrVmin_uV_clu, vrVmax_uV_clu, vrVpp_uV_clu] = clu_amp_(S_clu.tmrWav_raw_clu, S_clu.viSite_clu);
    mrVmin_clu = squeeze_(min(S_clu.tmrWav_clu));
    mrVmax_clu = squeeze_(max(S_clu.tmrWav_clu));
    mrVmin_uv_clu = squeeze_(min(S_clu.tmrWav_raw_clu));
    mrVmax_uv_clu = squeeze_(max(S_clu.tmrWav_raw_clu));
    % tmrWav_clu = S_clu.tmrWav_clu;
    % mrVmin_clu = shiftdim(min(tmrWav_clu,[],1));
    % mrVmax_clu = shiftdim(max(tmrWav_clu,[],1));
    vrVpp_clu = mr2vr_sub2ind_(mrVmax_clu-mrVmin_clu, S_clu.viSite_clu, 1:S_clu.nClu);
    vrVmin_clu = mr2vr_sub2ind_(mrVmin_clu, S_clu.viSite_clu, 1:S_clu.nClu);
    vrVpp_uv_clu = mr2vr_sub2ind_(mrVmax_uv_clu-mrVmin_uv_clu, S_clu.viSite_clu, 1:S_clu.nClu);
    vrVmin_uv_clu = mr2vr_sub2ind_(mrVmin_uv_clu, S_clu.viSite_clu, 1:S_clu.nClu);
    % [vrVpp_clu, ~] = max(mrVmax_clu - mrVmin_clu,[],1);
    % [vrVmin_clu, viSite_min_clu] = min(mrVmin_clu,[],1);
    % if ~isfield(S_clu, 'viSite_clu')
    %     viSite_clu = viSite_min_clu;
    % else
    %     viSite_clu = S_clu.viSite_clu;
    % end
    % vrVpp_clu=vrVpp_clu(:);  viSite_clu=viSite_clu(:);
    % [vrVpp_clu, viSite_clu, vrVmin_clu] = multifun_(@(x)x(:), vrVpp_clu, viSite_clu, vrVmin_clu);

    try
        S0 = get(0, 'UserData');
        if isempty(S0)
            S0 = load0_(subsFileExt_(P.vcFile_prm, '_jrc.mat'));
            set(0, 'UserData', S0);
        end
        vrVrms_site = bit2uV_(single(S0.vrThresh_site(:)) / S0.P.qqFactor, P);
        %     vrSnr_clu = vrVpp_clu ./ vrVrms_site(viSite_clu);
        vrSnr_clu = abs(vrVmin_clu) ./ vrVrms_site(S_clu.viSite_clu);
        vnSite_clu = sum(bsxfun(@lt, mrVmin_clu, -vrVrms_site * S0.P.qqFactor),1)';
    catch ME
        [vrVrms_site, vrSnr_clu, vnSite_clu] = deal([]);
        warning(ME.identifier, 'Could not compute vrVrms_site, vrSnr_clu, vnSite_clu: %s', ME.message);
    end
    [vrIsoDist_clu, vrLRatio_clu, vrIsiRatio_clu] = S_clu_quality2_(S_clu, P, viClu_update);

    S_clu = struct_add_(S_clu, vrVpp_clu, vrSnr_clu, vrVrms_site, vnSite_clu, ...
    vrIsoDist_clu, vrLRatio_clu, vrIsiRatio_clu, vrVpp_uv_clu, vrVmin_uv_clu);
    fprintf('\ttook %0.1fs\n', toc(t1));
end %func
