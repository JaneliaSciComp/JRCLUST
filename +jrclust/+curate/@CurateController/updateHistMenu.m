function updateHistMenu(obj)
    %UPDATEHISTMENU Update the history menu
    if obj.hasMenu('HistMenu')
        hClust = obj.hClust;
        hHistMenu = obj.hMenus('HistMenu');

        % show at most 5 entries in the history menu
        startAt = min(5, hClust.nEdits);

        delete(hHistMenu.Children);
        for i = startAt:-1:1
            label = obj.hClust.history.message{i};
            uimenu(hHistMenu, 'Label', label, ...
                   'Callback', @(hO, hE) obj.restoreHistory(i-1), ...
                   'Checked', jrclust.utils.ifEq(i == 1, 'on', 'off'), ...
                   'Enable', 'on');
        end
    end
end

