function probe(obj, probeFile)
    %PROBE Plot a probe layout
    if nargin > 1
        [~, ~, ext] = fileparts(probeFile);
        if isempty(ext) % a convenience for a forgetful mind
            probeFile = [probeFile '.prb'];
        end

        probeFile_ = jrclust.utils.absPath(probeFile, fullfile(jrclust.utils.basedir(), 'probes'));
        if isempty(probeFile_)
            obj.errMsg = sprintf('Could not find probe file: %s', probeFile);
            obj.isError = 1;
            return;
        end

        showProbe(jrclust.Config(struct('probe_file', probeFile)));
    elseif isempty(obj.hCfg) 
        obj.errMsg = 'Specify a probe file or config file';
        obj.isError = 1;
        return;
    else
        showProbe(obj.hCfg);
    end
end

%% LOCAL FUNCTIONS
function showProbe(probeData)
    %DOPLOTPROBE Plot a figure representing a probe
    hFigProbe = jrclust.views.Figure('FigProbe', [0 0 .5 1], 'Probe', 0, 0);
    vrX = 0.5*[-1 -1 1 1] * probeData.probePad(2);
    vrY = 0.5*[-1 1 1 -1] * probeData.probePad(1);

    uShanks = unique(probeData.shankMap);
    nShanks = numel(uShanks);
    [nRows, nCols] = min(ceil(nShanks ./ (1:4)), [], 2);

    hFigProbe.addSubplot('hShanks', nRows, nCols);

    for iShank = 1:nShanks
        shank = uShanks(iShank);
        shankMask = (probeData.shankMap == shank);

        XData = bsxfun(@plus, probeData.siteLoc(shankMask, 1)', vrX(:));
        YData = bsxfun(@plus, probeData.siteLoc(shankMask, 2)', vrY(:));
        nSites = sum(shankMask);

        % plot sites
        hFigProbe.subplotApply('hShanks', iShank, @patch, XData, YData, 'w', 'EdgeColor', 'k');

        % label sites
        if ~isempty(probeData.siteMap(shankMask))
            siteLabels = arrayfun(@(i) sprintf('%d/%d', i, probeData.siteMap(i)), find(shankMask), 'UniformOutput', 0);
        else
            siteLabels = arrayfun(@(i) sprintf('%d', i), 1:nSites, 'UniformOutput', 0);
        end

        hFigProbe.subplotApply('hShanks', iShank, @text, probeData.siteLoc(shankMask, 1), ...
                               probeData.siteLoc(shankMask, 2), ...
                               siteLabels, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

        spTitle = 'Site # / Chan #';
        if nShanks > 1
            spTitle = sprintf('Shank %d: %s', shank, spTitle);
        end

        hFigProbe.subplotApply('hShanks', iShank, @axis, [min(XData(:)), max(XData(:)), min(YData(:)), max(YData(:))]);
        hFigProbe.subplotApply('hShanks', iShank, @title, spTitle);
        hFigProbe.subplotApply('hShanks', iShank, @xlabel, 'X Position (\mum)');
        hFigProbe.subplotApply('hShanks', iShank, @ylabel, 'Y Position (\mum)');
        hFigProbe.subplotApply('hShanks', iShank, @axis, 'equal');
    end

    hFigProbe.setMouseable();
end