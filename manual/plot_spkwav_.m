%--------------------------------------------------------------------------
function S_fig = plot_spkwav_(S_fig, S0)
    % fPlot_raw = 0;
    if nargin<2, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    [P, spikeSites, S_clu] = deal(S0.P, S0.spikeSites, S0.S_clu);
    spikeWaveforms = get_spkwav_(P);

    [cvrX, cvrY, cviSite] = deal(cell(S_clu.nClusters, 1));
    vnSpk = zeros(S_clu.nClusters, 1);
    miSites_clu = P.miSites(:, S_clu.clusterSites);
    if isfield(S_fig, 'maxAmp')
        maxAmp = S_fig.maxAmp;
    else
        maxAmp = P.maxAmp;
    end
    for iCluster = 1:S_clu.nClusters
        try
            spikesToShow = randomSelect_(S_clu_viSpk_(S_clu, iCluster, spikeSites), P.nSpk_show);
            if P.fWav_raw_show
                trWav1 = raw2uV_(spikeWaveforms(:,:,spikesToShow), P);
                trWav1 = fft_lowpass_(trWav1, getOr(P, 'fc_spkwav_show', []), P.sampleRateHz);
            else
                trWav1 = tnWav2uV_(spikeWaveforms(:,:,spikesToShow), P);
            end
            sitesToShow = miSites_clu(:, iCluster);
            [cvrY{iCluster}, cvrX{iCluster}] = tr2plot_(trWav1, iCluster, sitesToShow, maxAmp, P);
            cviSite{iCluster} = sitesToShow;
            vnSpk(iCluster) = size(trWav1, 3); %subsample
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
