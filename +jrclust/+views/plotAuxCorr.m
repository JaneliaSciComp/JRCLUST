function corrData = plotAuxCorr(hClust, selected)
    %PLOTAUXCORR Plot aux channel correlation with firing rates figure
    if nargin < 2
        selected = [];
    end
    hCfg = hClust.hCfg;

    % load aux channel data
    [auxSamples, auxTimes] = loadAuxChannel(hCfg);
    if isempty(auxSamples)
        corrData = [];
        jrclust.utils.qMsgBox('Aux input not found');
        return;
    end

    % compute firing rates and correlate with the aux channel
    firingRates = hClust.getFiringRates([], numel(auxSamples));
    auxChanCorr = arrayfun(@(i) corr(auxSamples, firingRates(:, i), 'type', 'Pearson'), 1:size(firingRates, 2));

    [~, argsort] = sort(auxChanCorr, 'descend');

    nClustersShow = min(hCfg.nClustersShowAux, numel(auxChanCorr));
    auxLabel = hCfg.getOr('auxLabel', 'aux');
    nSubsamplesAux = hCfg.getOr('nSubsamplesAux', 100);

    if ~isempty(selected)
        nClustersShow = 1;
        argsort = selected;
    end

    hFigAux = jrclust.views.Figure('FigAux', [.5 0 .5 1], hCfg.sessionName, 1, 1);
    hTabGroup = hFigAux.figApply(@uitabgroup);

    for iiCluster = 1:nClustersShow
        iCluster = argsort(iiCluster);
        hTab = uitab(hTabGroup, 'Title', sprintf('Cluster %d', iCluster), 'BackgroundColor', 'w');
        axes('Parent', hTab);
        subplot(2, 1, 1);

        hAx = plotyy(auxTimes, firingRates(:, iCluster), auxTimes, auxSamples);

        xlabel('Time (s)');
        ylabel(hAx(1),'Firing Rate (Hz)');
        ylabel(hAx(2), auxLabel);

        iSite = hClust.clusterSites(iCluster);
        iTitle = sprintf('Cluster %d (Site %d, Chan %d): Corr=%0.3f', iCluster, iSite, hCfg.siteMap(iSite), auxChanCorr(iCluster));
        title(iTitle);
        set(hAx, 'XLim', auxTimes([1,end]));
        grid on;

        subplot(2, 1, 2);
        plot(auxSamples(1:nSubsamplesAux:end), firingRates(1:nSubsamplesAux:end,iCluster), 'k.');
        xlabel(auxLabel);
        ylabel('Firing Rate (Hz)');
        grid on;
    end

    corrData = struct('firingRates', firingRates, ...
                      'auxSamples', auxSamples, ...
                      'auxChanCorr', auxChanCorr, ...
                      'auxTimes', auxTimes);
end

%% LOCAL FUNCTIONS
function [auxSamples, auxTimes] = loadAuxChannel(hCfg)
    %LOADAUXCHANNEL Load the aux channel
    [auxSamples, auxTimes] = deal([]);
    if numel(hCfg.rawRecordings) > 1
        jrclust.utils.qMsgBox('Multi-file mode is currently not supported');
        return;
    end

    % try to guess auxFile
    if isempty(hCfg.auxFile)
        [~, ~, ext] = fileparts(hCfg.rawRecordings{1});

        if strcmpi(ext, '.ns5')
            try
                hCfg.auxFile = jrclust.utils.subsExt(hCfg.rawRecordings{1}, '.ns2');
            catch
                return;
            end
        elseif ismember(lower(ext), {'.bin', '.dat'})
            try
                hCfg.auxFile = hCfg.rawRecordings{1};
            catch
                return;
            end
        else
            return;
        end
    end

    [~, ~, auxExt] = fileparts(hCfg.auxFile);
    switch lower(auxExt)
%         case '.ns2'
%             auxChan = hCfg.getOr('auxChan', 1);
%             [mnWav_aux, hFile_aux, auxData] = load_nsx_(hCfg.auxFile);
%             scale_aux = hFile_aux.Entity(auxChan).Scale * hCfg.auxScale;
%             vrWav_aux = single(mnWav_aux(auxChan,:)') * scale_aux;
%             auxRate = auxData.sRateHz;
        case '.mat'
            auxData = load(hCfg.auxFile);
            auxDataFields = fieldnames(auxData);
            auxSamples = auxData.(auxDataFields{1});
            auxRate = hCfg.getOr('auxRate', hCfg.sampleRate);
        case {'.dat', '.bin'}
            if isempty(hCfg.auxChan)
                return;
            end

            hRec = jrclust.detect.newRecording(hCfg.auxFile, hCfg);
            auxSamples = single(hRec.readRawROI(hCfg.auxChan, 1:hRec.nSamples))*hCfg.bitScaling*hCfg.auxScale;
            auxRate = hCfg.getOr('auxRate', hCfg.sampleRate);
        otherwise
            jrclust.utils.qMsgBox(sprintf('hCfg.auxFile: unsupported file format: %s\n', auxExt));
        return;
    end % switch

    if nargout >= 2
        auxTimes = single(1:numel(auxSamples))'/auxRate;
    end
end
