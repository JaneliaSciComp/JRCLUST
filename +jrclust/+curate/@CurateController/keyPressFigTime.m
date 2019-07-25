function keyPressFigTime(obj, ~, hEvent)
    %KEYPRESSFIGTIME Handle callbacks for keys pressed in time view
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    hFigTime = obj.hFigs('FigTime');
    factor = 4^double(jrclust.utils.keyMod(hEvent, 'shift')); % 1 or 4

    switch hEvent.Key
        case 'leftarrow' % go down one channel
            obj.currentSite = max(obj.currentSite - factor, 1);
            obj.updateFigTime(0);

        case 'rightarrow' % go up one channel
            obj.currentSite = min(obj.currentSite + factor, obj.hCfg.nSites);
            obj.updateFigTime(0);

        case 'uparrow'
            jrclust.views.rescaleFigTime(hFigTime, sqrt(2)^-factor);

        case 'downarrow' % change amp
            jrclust.views.rescaleFigTime(hFigTime, sqrt(2)^factor);

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
            obj.splitPoly(hFigTime, jrclust.utils.keyMod(hEvent, 'shift'));

    end % switch
end