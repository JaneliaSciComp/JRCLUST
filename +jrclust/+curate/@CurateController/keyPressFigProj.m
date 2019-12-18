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
            if obj.projSites(1) - factor > 0
                obj.projSites = obj.projSites - factor;
            else
                obj.projSites = obj.projSites - obj.projSites(1) + 1;
            end
            obj.updateFigProj(0);

        case 'rightarrow' % go up one channel
            if obj.projSites(end) + factor <= obj.hCfg.nSites
                obj.projSites = obj.projSites + factor;
            else
                obj.projSites = obj.hCfg.nSites - (obj.hCfg.nSitesFigProj-1:-1:0);
            end
            obj.updateFigProj(0);

        case 'b' % background spikes
            hFigProj.toggleVisible('background');
            
        case {'d', 'backspace', 'delete'} % delete
            hFigProj.wait(1);
            obj.deleteClusters();
            hFigProj.wait(0);

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
            if ismember(obj.hCfg.dispFeature, {'pca', 'ppca', 'gpca', 'template'})
                % [1, 2] => [1, 3] => [2, 3] => [1, 2] => ...
                obj.hCfg.pcPair = sort(mod(obj.hCfg.pcPair + 1, 3) + 1);
                obj.updateFigProj(0);
            end

        case 'r' %reset view
            obj.updateFigProj(1);

        case 's' %split
            obj.splitPoly(hFigProj, jrclust.utils.keyMod(hEvent, 'shift'));
        
    end % switch
end