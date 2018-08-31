%--------------------------------------------------------------------------
function S_fig = plot_tnWav_clu_(S_fig, P)
    % Substituting plot_spk_
    S0 = get(0, 'UserData');
    S_clu = S0.S_clu;
    if ~isfield(P, 'LineWidth'), P.LineWidth=1; end
    trWav_clu = ifeq_(P.fWav_raw_show, S_clu.tmrWav_raw_clu, S_clu.tmrWav_clu);
    [nSamples, nSites, nClu] = size(trWav_clu);
    nChans_show = size(P.miSites, 1);
    miSites_clu = P.miSites(:, S_clu.clusterSites);
    % nSites = numel(P.chanMap);

    % determine x
    x_offset = P.spkLim(2) / (diff(P.spkLim)+1); %same for raw and filt
    vrX = (1:nSamples*nClu)/nSamples + x_offset;
    vrX(1:nSamples:end) = nan;
    vrX(nSamples:nSamples:end) = nan;
    trWav_clu = trWav_clu / S_fig.maxAmp;

    % nChans_show = size(P.miSites,1);
    mrX = repmat(vrX(:), [1, nChans_show]);
    mrX = reshape(mrX, [nSamples, nClu, nChans_show]);
    mrX = reshape(permute(mrX, [1 3 2]), [nSamples*nChans_show, nClu]);

    mrY = zeros(nSamples * nChans_show, nClu, 'single');
    for iClu=1:nClu
        viSites1 = miSites_clu(:,iClu);
        mrY1 = trWav_clu(:,viSites1,iClu);
        mrY1 = bsxfun(@plus, mrY1, single(viSites1'));
        mrY(:,iClu) = mrY1(:);
    end

    % if ~isfield(P, 'LineStyle') || isempty(P.LineStyle)
    if isfield(S_fig, 'vhPlot')
        plot_update_(S_fig.vhPlot, mrX, mrY);
    else
        S_fig.vhPlot = plot_group_(S_fig.hAx, mrX, mrY, 'LineWidth', P.LineWidth);
    end
    % else
    %     S_fig.vhPlot = plot_group_(S_fig.hAx, mrX, mrY, P.LineStyle, 'LineWidth', P.LineWidth);
    % end
    set(S_fig.hAx, 'YTick', 1:nSites, 'XTick', 1:nClu);
end % function
