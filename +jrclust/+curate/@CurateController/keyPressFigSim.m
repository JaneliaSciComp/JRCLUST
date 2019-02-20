function keyPressFigSim(obj, ~, hEvent)
    %KEYPRESSFIGSIM Handle callbacks for keys pressed in sim view
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    hFigSim = obj.hFigs('FigSim');

    switch hEvent.Key
        case {'d', 'backspace', 'delete'} % delete
            hFigSim.wait(1);
            obj.deleteClusters();
            hFigSim.wait(0);

        case 'k' % kilosort template similarity view
            if isa(obj.hClust, 'jrclust.sort.TemplateClustering') && strcmp(hFigSim.figData.figView, 'waveform')
                hFigSim.figData.figView = 'template';
                obj.updateFigSim();
            end

        case 'm' % merge
            hFigSim.wait(1);
            obj.mergeSelected();
            hFigSim.wait(0);

        case 's' % split
            hFigSim.wait(1);
            obj.autoSplit(1);
            hFigSim.wait(0);

        case 'w' % waveform-based sim view
            if isa(obj.hClust, 'jrclust.sort.TemplateClustering') && strcmp(hFigSim.figData.figView, 'template')
                hFigSim.figData.figView = 'waveform';
                obj.updateFigSim();
            end
    end % switch
end