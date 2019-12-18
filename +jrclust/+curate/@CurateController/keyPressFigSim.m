function keyPressFigSim(obj, ~, hEvent)
    %KEYPRESSFIGSIM Handle callbacks for keys pressed in sim view
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    hFigSim = obj.hFigs('FigSim');

    switch hEvent.Key
        case 'uparrow'
            selected = obj.selected;
            if numel(selected) == 1
                selected = [selected selected];
            end
            selected(2) = obj.showSubset(min(obj.nShown, obj.unitIndex(selected(2)) + 1));
            obj.updateSelect(selected);

        case 'downarrow'
            selected = obj.selected;
            if numel(selected) == 1
                selected = [selected selected];
            end
            selected(2) = obj.showSubset(max(1, obj.unitIndex(selected(2)) - 1));
            obj.updateSelect(selected);

        case 'rightarrow'
            selected = obj.selected;
            if numel(selected) == 1
                selected = [selected selected];
            end
            selected(1) = obj.showSubset(min(obj.nShown, obj.unitIndex(selected(1)) + 1));
            obj.updateSelect(selected);

        case 'leftarrow'
            selected = obj.selected;
            if numel(selected) == 1
                selected = [selected selected];
            end
            selected(1) = obj.showSubset(max(1, obj.unitIndex(selected(1)) - 1));
            obj.updateSelect(selected);
            
        case {'d', 'backspace', 'delete'} % delete
            hFigSim.wait(1);
            obj.deleteClusters();
            hFigSim.wait(0);

        case 'h' % help
            jrclust.utils.qMsgBox(hFigSim.figData.helpText, 1);

        case 'k' % kilosort template similarity view
            if isa(obj.hClust, 'jrclust.sort.TemplateClustering') && strcmp(hFigSim.figData.figView, 'waveform')
                hFigSim.figData.figView = 'template';
                obj.updateFigSim();
            end

        case 'm' % merge
            hFigSim.wait(1);
            obj.mergeSelected();
            hFigSim.wait(0);
            
        case 'r' %reset view
            obj.updateCursorFigSim();

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