function updateFigRD(obj)
    %UPDATEFIGRD Update the rho-delta plot
    if ~obj.hasFig('FigRD')
        return;
    end

    jrclust.views.plotFigRD(obj.hFigs('FigRD'), obj.hClust, obj.hCfg);
end
