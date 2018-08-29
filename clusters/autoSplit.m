%--------------------------------------------------------------------------
function autoSplit(fMulti, S0)
    % Auto-split feature that calls Hidehiko Inagaki's code
    % 20160426
    if nargin < 1
        fMulti = 0;
    end
    if nargin < 2
        S0 = [];
    end
    if isempty(S0)
        S0 = get(0, 'UserData');
    end
    P = S0.P;
    S_clu = S0.S_clu;

    if ~isempty(S0.secondarySelectedCluster)
        msgbox_('Select one cluster');
        return;
    end
    if S_clu.nSpikesPerCluster(S0.primarySelectedCluster) < 2
        msgbox_('At least two spikes required for splitting');
        return;
    end

    hMsg = msgbox_('Splitting... (this closes automatically)');
    clusterToSplit = S0.primarySelectedCluster;
    clusterSite = S_clu.clusterSites(clusterToSplit);

    if fMulti
        clusterSites = P.miSites(1:end-P.nSites_ref, clusterSite);
    else
        clusterSites = clusterSite;
    end

    tmp = getSpikeWaveformsSites(S_clu.spikesByCluster{clusterToSplit}, clusterSites, S0);
    mrSpkWav1 = tnWav2uV_(tmp, P);
    mrSpkWav1 = reshape(mrSpkWav1, [], size(mrSpkWav1,3));

    [~, temp_S_fig] = getCachedFig('FigTime'); % TW gets the variables in the "time" figure -- needed for getting the highlighted site
    siteToUse = temp_S_fig.iSite; % TW use the highlighted site from "time" figure
    % siteToUse = clusterSite; % TW use dominant site for the cluster

    % TW calculate amplitudes on the fly
    mrWav_spk1 = squeeze_(tnWav2uV_(getSpikeWaveformsSites(S_clu.spikesByCluster{clusterToSplit}, siteToUse, S0), P));
    mrFet1 = max(mrWav_spk1) - min(mrWav_spk1);

    % [vlSpkIn, mrFet_split, vhAx] = auto_split_wav_(mrSpkWav1, [S0.spikeTimes(S_clu.spikesByCluster{clusterToSplit}) S0.vrAmp_spk(S_clu.spikesByCluster{clusterToSplit})], 2); % MNE
    [vlSpkIn, mrFet_split, vhAx] = auto_split_wav_(mrSpkWav1, [S0.spikeTimes(S_clu.spikesByCluster{clusterToSplit}) mrFet1'], 2); %TW

    hPoly = [];
    hFigTemp = gcf;
    try, drawnow; close(hMsg); catch; end
    while 1
        vcAns = userDialog('Split?', 'confirmation', 'Yes', 'No', 'Manual', 'Yes');
        switch lower(vcAns)
            case 'yes'
                close(hFigTemp);
                break;

            case {'no', ''}
                close(hFigTemp);
                return;

            case 'manual'
                vcAns = userDialog('Select projection', '', 'PC1 vs PC2', 'PC3 vs PC2', 'PC1 vs PC3', 'PC1 vs PC2');
                switch vcAns
                    case 'PC1 vs PC2'
                        [hAx_, iAx1, iAx2] = deal(vhAx(1), 1, 2);

                    case 'PC3 vs PC2'
                        [hAx_, iAx1, iAx2] = deal(vhAx(2), 3, 2);

                    case 'PC1 vs PC3',
                        [hAx_, iAx1, iAx2] = deal(vhAx(3), 1, 3);

                    otherwise
                        close(hFigTemp); return;
                end

                axes(hAx_);
                cla(hAx_);

                [vrX1, vrY1] = deal(mrFet_split(:, iAx1), mrFet_split(:, iAx2));
                plot(hAx_, vrX1, vrY1, 'k.');
                tryClose(hPoly);
                hPoly = impoly_();
                mrPolyPos = getPosition(hPoly);
                vlSpkIn = inpolygon(vrX1, vrY1, mrPolyPos(:,1), mrPolyPos(:,2));
                plot(vrX1(vlSpkIn), vrY1(vlSpkIn), 'b.', vrX1(~vlSpkIn), vrY1(~vlSpkIn), 'r.');
        end % switch
    end

    splitCluster(clusterToSplit, vlSpkIn);
end % function
