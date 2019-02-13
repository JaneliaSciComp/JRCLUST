function keyPressFigTime(obj, ~, hEvent)
    %KEYPRESSFIGTIME Handle callbacks for keys pressed in time view
    if obj.isWorking
        return;
    end

    hFigTime = obj.hFigs('FigTime');
    factor = 4^double(keyMod(hEvent, 'shift')); % 1 or 4

    switch hEvent.Key
        case 'leftarrow' % go down one channel
            obj.currentSite = max(obj.currentSite - factor, 1);
            obj.updateFigTime(0);

        case 'rightarrow' % go up one channel
%             if ~isVisible_(S_fig.hAx)
%                 jrclust.utils.qMsgBox('Channel switching is disabled in the position view'); return;
%             end
            obj.currentSite = min(obj.currentSite + factor, obj.hCfg.nSites);
            obj.updateFigTime(0);

        case 'uparrow'
            rescaleFigTime(hFigTime, sqrt(2)^-factor);

        case 'downarrow' % change amp
            rescaleFigTime(hFigTime, sqrt(2)^factor);

        case 'b' % toggle background spikes
            hFigTime.figData.doPlotBG = hFigTime.toggleVisible('background');

        case 'f' % toggle feature display
            if strcmp(obj.hCfg.dispFeature, 'vpp')
                obj.updateProjection(obj.hCfg.clusterFeature);
            else
                obj.updateProjection('vpp');
            end

        case 'h' % help
            jrclust.utils.qMsgBox(hFigTime.figData.helpText, 1);

        case 'm' % merge
            hFigTime.wait(1);
            obj.mergeSelected();
            hFigTime.wait(0);

        case 'r' % reset view
            obj.updateFigTime(1);

        case 's' % split
            if numel(obj.selected) == 1
                iCluster = obj.selected(1);

                hFigTime.addPlot('hPoly', @impoly)
                try
                    polyPos = hFigTime.plotApply('hPoly', @getPosition);
                catch ME
                    return;
                end

                XData = hFigTime.plotApply('foreground', @get, 'XData');
                YData = hFigTime.plotApply('foreground', @get, 'YData');

                retained = inpolygon(XData, YData, polyPos(:,1), polyPos(:,2));
                hFigTime.addPlot('hSplit', @line, XData(retained), YData(retained), ...
                                 'Color', obj.hCfg.colorMap(3, :), ...
                                 'Marker', '.', 'LineStyle', 'none');

                dlgAns = questdlg('Split?', 'Confirmation', 'No');

                hFigTime.rmPlot('hPoly');
                hFigTime.rmPlot('hSplit');

                if strcmp(dlgAns, 'Yes')
                    obj.splitCluster(iCluster, retained);
                end
            end
    end % switch
end