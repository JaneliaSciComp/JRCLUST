function plotAllFigures(obj)
    %PLOTALLFIGURES Plot all figures
    if isempty(obj.hFigs)
        obj.spawnFigures();
    elseif ~all(obj.figApply(@(hFig) hFig.isReady)) % clean up from an aborted session
        obj.closeFigures();
        obj.spawnFigures();
    end

    % plot sim score figure
    if obj.hasFig('FigSim')
        % set key and mouse handles
        hFigSim = obj.hFigs('FigSim');
        obj.updateFigSim();
        hFigSim.hFunKey = @obj.keyPressFigSim;
        hFigSim.setMouseable(@obj.mouseClickFigSim);
    end

    % plot feature projection
    if obj.hasFig('FigProj')
        % set key and mouse handles
        hFigProj = obj.hFigs('FigProj');
        hFigProj.hFunKey = @obj.keyPressFigProj;
        hFigProj.setMouseable(); % no special mouse function
    end

    % plot amplitude vs. time
    if obj.hasFig('FigTime')
        % set key and mouse handles
        hFigTime = obj.hFigs('FigTime');
        hFigTime.hFunKey = @obj.keyPressFigTime;
    end

    % plot main waveform view
    if obj.hasFig('FigWav')
        % set key and mouse handles
        hFigWav = jrclust.views.plotFigWav(obj.hFigs('FigWav'), obj.hClust, obj.hCfg, obj.maxAmp);
        setFigWavXTicks(hFigWav, obj.hClust, 1); % show cluster counts by default

        hFigWav.hFunKey = @obj.keyPressFigWav;
        hFigWav.setMouseable(@obj.mouseClickFigWav);

        % make this guy the key log
        hFigWav.figApply(@set, 'CloseRequestFcn', @obj.killFigWav);
        obj.addMenu(hFigWav);
    end

    % plot rho-delta figure
    obj.updateFigRD();

    % select first cluster (also plots other figures)
    obj.updateSelect(1);

    % zoom in on first cluster
    obj.keyPressFigWav([], struct('Key', 'z')); % zoom in
end