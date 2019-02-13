function closeFigures(obj)
    %CLOSEFIGURES Close all open figures
    if obj.hasFig('FigWav') && obj.hFigs('FigWav').isReady
        hFigWav = obj.hFigs('FigWav');
        remove(obj.hFigs, 'FigWav'); % to prevent infinite recursion!
        hFigWav.close(); % calls killFigWav
    end

    try
        obj.figApply(@(hFig) hFig.close());
    catch ME
        warning('Could not close figures: %s', ME.message);
    end

    obj.hFigs = containers.Map();
end