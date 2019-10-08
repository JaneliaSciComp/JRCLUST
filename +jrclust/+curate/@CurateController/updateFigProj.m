function updateFigProj(obj, doAutoscale)
    %UPDATEFIGPROJ
    if ~obj.hasFig('FigProj')
        return;
    end

    hFigProj = obj.hFigs('FigProj');

    if isempty(hFigProj.figData) || ~isfield(hFigProj.figData, 'boundScale')
        boundScale = obj.maxAmp;
    else
        boundScale = hFigProj.figData.boundScale;
    end

    jrclust.views.plotFigProj(hFigProj, obj.hClust, obj.projSites, obj.selected, boundScale, obj.channel_idx, doAutoscale);

    hFigProj.setMouseable(); % no special mouse function
end
