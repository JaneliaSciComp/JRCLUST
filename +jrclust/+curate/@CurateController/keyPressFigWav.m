function keyPressFigWav(obj, ~, hEvent)
    %KEYPRESSFIGWAV Handle callbacks for keys pressed in main view
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    hFigWav = obj.hFigs('FigWav');
    factor = 4^double(jrclust.utils.keyMod(hEvent, 'shift')); % 1 or 4
    nSites = obj.hCfg.nSites;

    switch hEvent.Key
        case 'uparrow'
            obj.maxAmp = jrclust.views.rescaleFigWav(hFigWav, obj.hClust, obj.hCfg, obj.maxAmp, sqrt(2)^-factor);
            obj.updateCursorFigWav();

        case 'downarrow'
            obj.maxAmp = jrclust.views.rescaleFigWav(hFigWav, obj.hClust, obj.hCfg, obj.maxAmp, sqrt(2)^factor);
            obj.updateCursorFigWav();

        case 'leftarrow' % select previous cluster
            if jrclust.utils.keyMod(hEvent, 'shift')
                selected_ = [obj.selected(1), max(obj.selected(end)-1, 1)];
            else
                selected_ = max(obj.selected(1)-1, 1);
            end
            obj.updateSelect(selected_);

        case 'rightarrow' % select next cluster
            if jrclust.utils.keyMod(hEvent, 'shift')
                selected_ = [obj.selected(1), min(obj.selected(end)+1, obj.hClust.nClusters)];
            else
                selected_ = min(obj.selected(1)+1, obj.hClust.nClusters);
            end
            obj.updateSelect(selected_);

        case 'home' % select first cluster
            obj.updateSelect(1);
            obj.keyPressFigWav([], struct('Key', 'z')); % zoom in

        case 'end' % select last cluster
            obj.updateSelect(obj.hClust.nClusters);
            obj.keyPressFigWav([], struct('Key', 'z')); % zoom in

        case {'d', 'backspace', 'delete'}
            hFigWav.wait(1);
            obj.deleteClusters();
            hFigWav.wait(0);

        case 'm' % merge clusters
            hFigWav.wait(1);
            obj.mergeSelected();
            hFigWav.wait(0);

        case 'n' % toggle spike count in clusters
            obj.hCfg.showSpikeCount = ~obj.hCfg.showSpikeCount;
            setFigWavXTicks(hFigWav, obj.hClust, obj.hCfg.showSpikeCount);

        case 'space' % select most similar to currently selected
            waveformSim = obj.hClust.waveformSim;
            waveformSim(obj.selected(1), obj.selected(1)) = -inf;
            [~, nextBest] = max(waveformSim(:, obj.selected(1)));
            obj.updateSelect([obj.selected(1), nextBest]);

        case 'p' % PSTH plot
            if isempty(obj.hCfg.trialFile)
                jrclust.utils.qMsgBox('''trialFile'' not set. Reload .prm file after setting (under "File menu")');
                return;
            end

            obj.updateFigPSTH(1);

        case 's' % split
            hFigWav.wait(1);
            obj.autoSplit(1);
            hFigWav.wait(0);

        case 'r' %reset view
            hFigWav.wait(1);
            hFigWav.axApply('default', @axis, [0, obj.hClust.nClusters + 1, 0, obj.hCfg.nSites + 1]);
            hFigWav.wait(0);

        case 'w' % toggle individual spike waveforms
            hFigWav.toggleVisible('hSpkAll');

        case 'z' % zoom
            if isempty(obj.selected)
                obj.updateSelect(1);
            else
                iCluster = obj.selected(1);
                iSite = obj.hClust.clusterSites(iCluster);
                hFigWav.setWindow(iCluster + [-1, 1]*6, iSite + [-1, 1]*(obj.hCfg.nSiteDir*2+1), [0 obj.hClust.nClusters+1], [0 nSites+1]);
            end

        case 'a'
            obj.updateFigWav();
            obj.updateSelect(obj.selected);

        case 'h'
            jrclust.utils.qMsgBox(hFigWav.figData.helpText, 1);

        case {'0', 'numpad0'}
            obj.annotateUnit('to_delete', 0); % TW

        case {'1', 'numpad1'}
            obj.annotateUnit('single', 0); % TW

        case {'2', 'numpad2'}
            obj.annotateUnit('multi', 0); % TW

        otherwise
            hFigWav.wait(0); %stop waiting
    end
end