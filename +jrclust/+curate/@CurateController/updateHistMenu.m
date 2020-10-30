function updateHistMenu(obj)
%UPDATEHISTMENU Add entries to the history menu reflecting the last 5 edit
%operations performed by the user.
if obj.hasMenu('HistMenu')
    hClust = obj.hClust;
    hHistMenu = obj.hMenus('HistMenu');

    % show at most 5 entries in the history menu
    startAt = max(0, hClust.nEdits - 4);

    delete(hHistMenu.Children);
    for i = startAt:hClust.nEdits
        if i == 0
            label = 'Initial commit';
        else
            label = obj.hClust.history.message{i};
        end

        uimenu(hHistMenu, 'Label', label, ...
               'Callback', @(hO, hE) obj.revertLast(hClust.nEdits - i), ...
               'Checked', jrclust.utils.ifEq(i == hClust.nEdits, 'on', 'off'), ...
               'Enable', 'on');
    end
end
end

