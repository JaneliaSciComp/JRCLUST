%--------------------------------------------------------------------------
function auto_split_(fMulti, S0)
    % Auto-split feature that calls Hidehiko Inagaki's code
    % 20160426
    if nargin<1, fMulti = 0; end
    if nargin<2, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    [P, S_clu] = deal(S0.P, S0.S_clu);

    if ~isempty(S0.iCluPaste), msgbox_('Select one cluster'); return; end
    if S_clu.vnSpk_clu(S0.iCluCopy)<2, msgbox_('At least two spikes required for splitting'); return; end

    hMsg = msgbox_('Splitting... (this closes automatically)');
    iClu1 = S0.iCluCopy;
    iSite1 = S_clu.viSite_clu(iClu1);
    if fMulti
        viSites1 = P.miSites(1:end-P.nSites_ref, iSite1);
    else
        viSites1 = iSite1;
    end
    % mrSpkWav1 = tnWav2uV_(tnWav_sites_(tnWav_spk, S_clu.cviSpk_clu{iClu1}, viSites1));
    mrSpkWav1 = tnWav2uV_(tnWav_spk_sites_(S_clu.cviSpk_clu{iClu1}, viSites1, S0), P);
    mrSpkWav1 = reshape(mrSpkWav1, [], size(mrSpkWav1,3));

    [~, temp_S_fig] = getCachedFig('FigTime'); % TW gets the variables in the "time" figure -- needed for getting the highlighted site
    site_to_use = temp_S_fig.iSite; % TW use the highlighted site from "time" figure
    % site_to_use = iSite1; % TW use dominant site for the cluster
    mrWav_spk1 = squeeze_(tnWav2uV_(tnWav_spk_sites_(S_clu.cviSpk_clu{iClu1}, site_to_use, S0), P)); % TW calculate amplitudes on the fly
    mrFet1 = max(mrWav_spk1)-min(mrWav_spk1); % TW calculate amplitudes on the fly

    % [vlSpkIn, mrFet_split, vhAx] = auto_split_wav_(mrSpkWav1, [S0.spikeTimes(S_clu.cviSpk_clu{iClu1}) S0.vrAmp_spk(S_clu.cviSpk_clu{iClu1})], 2); % MNE
    [vlSpkIn, mrFet_split, vhAx] = auto_split_wav_(mrSpkWav1, [S0.spikeTimes(S_clu.cviSpk_clu{iClu1}) mrFet1'], 2); %TW

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
            close(hFigTemp); return;
            case 'manual'
            vcAns = userDialog('Select projection', '', 'PC1 vs PC2', 'PC3 vs PC2', 'PC1 vs PC3', 'PC1 vs PC2');
            switch vcAns
                case 'PC1 vs PC2', [hAx_, iAx1, iAx2] = deal(vhAx(1), 1, 2);
                case 'PC3 vs PC2', [hAx_, iAx1, iAx2] = deal(vhAx(2), 3, 2);
                case 'PC1 vs PC3', [hAx_, iAx1, iAx2] = deal(vhAx(3), 1, 3);
                otherwise
                close(hFigTemp); return;
            end
            %             msgbox_(sprintf('Draw a polygon in PC%d vs PC%d', iAx1, iAx2), 1);
            axes(hAx_); cla(hAx_);
            [vrX1, vrY1] = deal(mrFet_split(:,iAx1), mrFet_split(:,iAx2));
            plot(hAx_, vrX1, vrY1, 'k.');
            tryClose(hPoly);
            hPoly = impoly_();
            mrPolyPos = getPosition(hPoly);
            vlSpkIn = inpolygon(vrX1, vrY1, mrPolyPos(:,1), mrPolyPos(:,2));
            plot(vrX1(vlSpkIn), vrY1(vlSpkIn), 'b.', vrX1(~vlSpkIn), vrY1(~vlSpkIn), 'r.');
        end %switch
    end
    split_clu_(iClu1, vlSpkIn);
end %func
