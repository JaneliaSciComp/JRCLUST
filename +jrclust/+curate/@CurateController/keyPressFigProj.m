function keyPressFigProj(obj, ~, hEvent)
    %KEYPRESSFIGPROJ Handle callbacks for keys pressed in feature view
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    hFigProj = obj.hFigs('FigProj');
    factor = 4^double(jrclust.utils.keyMod(hEvent, 'shift')); % 1 or 4

    switch lower(hEvent.Key)
        case 'uparrow'
            projScale = hFigProj.figData.boundScale*sqrt(2)^-factor;
            jrclust.views.rescaleFigProj(hFigProj, projScale, obj.hCfg);

        case 'downarrow'
            projScale = hFigProj.figData.boundScale*sqrt(2)^factor;
            jrclust.views.rescaleFigProj(hFigProj, projScale, obj.hCfg);

        case 'leftarrow' % go down one channel
            if obj.projSites(1) > 1
                obj.projSites = obj.projSites - 1;
                obj.updateFigProj(0);
            end

        case 'rightarrow' % go up one channel
            if obj.projSites(end) < obj.hCfg.nSites
                obj.projSites = obj.projSites + 1;
                obj.updateFigProj(0);
            end

        case 'b' % background spikes
            hFigProj.toggleVisible('background');

        case 'f' % toggle feature display
            if strcmp(obj.hCfg.dispFeature, 'vpp')
                if isa(obj.hClust, 'jrclust.sort.TemplateClustering')
                    obj.updateProjection('template');
                else
                    obj.updateProjection(obj.hCfg.clusterFeature);
                end
            else
                obj.updateProjection('vpp');
            end

        case 'h' % help
            jrclust.utils.qMsgBox(hFigProj.figData.helpText, 1);

        case 'm' % merge clusters
            hFigProj.wait(1);
            obj.mergeSelected();
            hFigProj.wait(0);

        case 'p' % toggle PCi v. PCj
            if ismember(obj.hCfg.dispFeature, {'pca', 'template'})
                % [1, 2] => [1, 3] => [2, 3] => [1, 2] => ...
                obj.hCfg.pcPair = sort(mod(obj.hCfg.pcPair + 1, 3) + 1);
                obj.updateFigProj(0);
            end

        case 'r' %reset view
            obj.updateFigProj(1);

        case 's' %split
            if numel(obj.selected) == 1
                iCluster = obj.selected(1);

                hFigProj.addPlot('hPoly', @impoly)
                polyPos = hFigProj.plotApply('hPoly', @getPosition);
                if isempty(polyPos)
                    return;
                end

                XData = hFigProj.plotApply('foreground', @get, 'XData');
                YData = hFigProj.plotApply('foreground', @get, 'YData');

                retained = inpolygon(XData, YData, polyPos(:,1), polyPos(:,2));

                % return here
                hFigProj.addPlot('hSplit', @line, XData(retained), YData(retained), ...
                                 'Color', obj.hCfg.colorMap(3, :), ...
                                 'Marker', '.', 'LineStyle', 'none');

                dlgAns = questdlg('Split?', 'Confirmation', 'No');

                hFigProj.rmPlot('hPoly');
                hFigProj.rmPlot('hSplit');

                if strcmp(dlgAns, 'Yes')
                    iSpikes = obj.hClust.spikesByCluster{iCluster};
                    retained = ismember(hFigProj.figData.dispFeatures.fgSpikes(retained), iSpikes);
                    obj.splitCluster(iCluster, retained);
                end
            end
        
    end % switch
end