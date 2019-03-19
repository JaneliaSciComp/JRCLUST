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

        case 'space' % select most similar to currently selected
            waveformSim = obj.hClust.waveformSim;
            waveformSim(obj.selected(1), obj.selected(1)) = -inf;
            [~, nextBest] = max(waveformSim(:, obj.selected(1)));
            obj.updateSelect([obj.selected(1), nextBest]);

        case {'0', 'numpad0'}
            obj.annotateUnit('to_delete', 0); % TW

        case {'1', 'numpad1'}
            obj.annotateUnit('single', 0); % TW

        case {'2', 'numpad2'}
            obj.annotateUnit('multi', 0); % TW

        case 'a'
            obj.updateFigWav();
            obj.updateSelect(obj.selected);

        case {'d', 'backspace', 'delete'}
            hFigWav.wait(1);
            obj.deleteClusters();
            hFigWav.wait(0);

        case 'g'
            dlgAns = jrclust.utils.inputdlgNum('Go to a cluster', '', 1);
            if ~isnan(dlgAns) && dlgAns > 0 && dlgAns <= obj.hClust.nClusters
                obj.updateSelect(dlgAns);
            end

        case 'h'
            jrclust.utils.qMsgBox(hFigWav.figData.helpText, 1);

        case 'm' % merge clusters
            hFigWav.wait(1);
            obj.mergeSelected();
            hFigWav.wait(0);

        case 'n' % toggle spike count in clusters
            obj.hCfg.showSpikeCount = ~obj.hCfg.showSpikeCount;
            setFigWavXTicks(hFigWav, obj.hClust, obj.hCfg.showSpikeCount);

        case 'p' % PSTH plot
            if isempty(obj.hCfg.trialFile)
                jrclust.utils.qMsgBox('''trialFile'' not set. Reload .prm file after setting (under "File menu")');
                return;
            end

            obj.updateFigPSTH(1);

        case 'q' % show and export quality scores for selected cluster
            iCluster = obj.selected(1);
            qScores = struct('unitID', iCluster, ...
                             'unitCenterSite', obj.hClust.clusterSites(iCluster), ...
                             'unitCount', obj.hClust.unitCount(iCluster), ...
                             'unitCentroid', obj.hClust.clusterCentroids(iCluster, :), ...
                             'unitPeak', obj.hClust.unitPeaksRaw(iCluster), ...
                             'unitVpp', obj.hClust.unitVppRaw(iCluster), ...
                             'unitIsoDist', obj.hClust.unitIsoDist(iCluster), ...
                             'unitLRatio', obj.hClust.unitLRatio(iCluster), ...
                             'unitISIRatio', obj.hClust.unitISIRatio(iCluster), ...
                             'unitNote', obj.hClust.clusterNotes{iCluster});

            qText = {sprintf('Quality metrics for unit %d', iCluster), ...
                     sprintf('Center site: %d (Peak site number which contains the most negative peak amplitude)', qScores.unitCenterSite), ...
                     sprintf('Unit count: %d (Number of spikes in unit)', qScores.unitCount), ...
                     sprintf('Centroid: %s (x (width) and y (depth, from tip) center-of-mass)', jrclust.utils.field2str(qScores.unitCentroid)), ...
                     sprintf('unitPeak: %0.3f (Min. voltage (uV) of the mean raw waveforms at the peak site (microvolts))', qScores.unitPeak), ...
                     sprintf('unitVpp: %0.3f (Peak-to-peak voltage (microvolts))', qScores.unitVpp), ...
                     sprintf('unitIsoDist: %0.3f (Isolation distance)', qScores.unitIsoDist), ...
                     sprintf('unitLRatio: %0.3f (L-ratio)', qScores.unitLRatio), ...
                     sprintf('unitISIRatio: %0.3f (ISI-ratio)', qScores.unitISIRatio), ...
                     sprintf('unitNote: ''%s'' (User comments)', qScores.unitNote)};

            if isprop(obj, 'unitSNR')
                qScores.unitSNR = obj.unitSNR(iCluster);
                qText{end+1} = sprintf('SNR: %0.3f (|Vp/Vrms|; Vp: negative peak amplitude of the peak site; Vrms: SD of the Gaussian noise (estimated from MAD))', qScores.unitSNR);
            end

            jrclust.utils.qMsgBox(qText);

            if jrclust.utils.keyMod(hEvent, 'shift')
                jrclust.utils.exportToWorkspace(qScores, 0);
            end

        case 'r' % reset view
            hFigWav.wait(1);
            hFigWav.axApply('default', @axis, [0, obj.hClust.nClusters + 1, 0, obj.hCfg.nSites + 1]);
            hFigWav.wait(0);

        case 's' % split
            hFigWav.wait(1);
            obj.autoSplit(1);
            hFigWav.wait(0);

        case 'w' % toggle individual spike waveforms
            hFigWav.toggleVisible('hSpkAll');

        case 'z' % zoom
            if isempty(obj.selected)
                obj.updateSelect(1);
            else
                iCluster = obj.selected(1);
                iSite = obj.hClust.clusterSites(iCluster);

                % do we have a second selected cluster?
                if numel(obj.selected) > 1
                    xRange = iCluster + [-1, 1]*max(abs(diff(obj.selected)) + 1, 6);
                else
                    xRange = iCluster + [-1, 1]*6;
                end
                hFigWav.setWindow(xRange, iSite + [-1, 1]*(obj.hCfg.nSiteDir*2+1), [0 obj.hClust.nClusters+1], [0 nSites+1]);
            end

        otherwise
            hFigWav.wait(0); %stop waiting
    end
end