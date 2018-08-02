%--------------------------------------------------------------------------
function S_fig = plot_spkwav_(S_fig, S0)
    % fPlot_raw = 0;
    if nargin<2, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    [P, spikeSites, S_clu] = deal(S0.P, S0.spikeSites, S0.S_clu);
    tnWav = get_spkwav_(P);

    [cvrX, cvrY, cviSite] = deal(cell(S_clu.nClusters, 1));
    vnSpk = zeros(S_clu.nClusters, 1);
    miSites_clu = P.miSites(:, S_clu.clusterSites);
    if isfield(S_fig, 'maxAmp')
        maxAmp = S_fig.maxAmp;
    else
        maxAmp = P.maxAmp;
    end
    for iClu = 1:S_clu.nClusters
        try
            viSpk_show = randomSelect_(S_clu_viSpk_(S_clu, iClu, spikeSites), P.nSpk_show);
            if P.fWav_raw_show
                trWav1 = raw2uV_(tnWav(:,:,viSpk_show), P);
                trWav1 = fft_lowpass_(trWav1, get_set_(P, 'fc_spkwav_show', []), P.sRateHz);
            else
                trWav1 = tnWav2uV_(tnWav(:,:,viSpk_show), P);
            end
            viSite_show = miSites_clu(:, iClu);
            [cvrY{iClu}, cvrX{iClu}] = tr2plot_(trWav1, iClu, viSite_show, maxAmp, P);
            cviSite{iClu} = viSite_show;
            vnSpk(iClu) = size(trWav1, 3); %subsample
        catch
            disperr_();
        end
    end
    S = makeStruct_(cvrY, cviSite, vnSpk);
    try
        set(S_fig.hSpkAll, 'XData', cell2mat_(cvrX), 'YData', cell2mat_(cvrY), 'UserData', S);
    catch
        S_fig.hSpkAll = plot(S_fig.hAx, cell2mat_(cvrX), cell2mat_(cvrY), 'Color', [.5 .5 .5], 'LineWidth', .5); %, P.LineStyle);
        set(S_fig.hSpkAll, 'UserData', S);
    end
end %func
