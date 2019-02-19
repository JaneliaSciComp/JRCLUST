function updateNoteMenu(obj)
    %UPDATENOTEMENU Update menu entry to indicate selected clusters
    if numel(obj.selected) > 1 && obj.hasMenu('InfoMenu')
        menuLabel = sprintf('Unit %d "%s" vs. Unit %d "%s"', obj.selected(1), ...
                            obj.hClust.clusterNotes{obj.selected(1)}, obj.selected(2), ...
                            obj.hClust.clusterNotes{obj.selected(2)});
        set(obj.hMenus('InfoMenu'), 'Label', menuLabel);
    elseif obj.hasMenu('InfoMenu')
        menuLabel = sprintf('Unit %d "%s"', obj.selected(1), obj.hClust.clusterNotes{obj.selected(1)});
        set(obj.hMenus('InfoMenu'), 'Label', menuLabel);
    end
end