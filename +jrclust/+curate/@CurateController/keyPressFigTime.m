function keyPressFigTime(obj, ~, hEvent)
    %KEYPRESSFIGTIME Handle callbacks for keys pressed in time view
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    hFigTime = obj.hFigs('FigTime');
    factor = 4^double(jrclust.utils.keyMod(hEvent, 'shift')); % 1 or 4

    switch hEvent.Key
        case 'uparrow'
            jrclust.views.rescaleFigTime(hFigTime, sqrt(2)^-factor);

        case 'downarrow' % change amp
            jrclust.views.rescaleFigTime(hFigTime, sqrt(2)^factor);

        case 'leftarrow' % go down one channel
            obj.currentSite = obj.spatial_idx(max(obj.channel_idx(obj.currentSite) - factor, 1));
            obj.updateFigTime(1);

        case 'rightarrow' % go up one channel
            obj.currentSite = obj.spatial_idx(min(obj.channel_idx(obj.currentSite) + factor, max(obj.channel_idx)));
            obj.updateFigTime(1);

        case 'b' % toggle background spikes
            hFigTime.figData.doPlotBG = hFigTime.toggleVisible('background');
            hFigTime.toggleVisible('background_hist');
            
        case {'d', 'backspace', 'delete'} % delete
            hFigTime.wait(1);
            obj.deleteClusters();
            hFigTime.wait(0);

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