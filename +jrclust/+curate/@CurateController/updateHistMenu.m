function updateHistMenu(obj)
    %UPDATEHISTMENU Update the history menu
    if obj.hasMenu('HistMenu')
        hClust = obj.hClust;
        hHistMenu = obj.hMenus('HistMenu');

        endAt = hClust.nEdits;
        startAt = max(0, endAt - 4); % show at most 5 (most recent) entries in hist menu

        delete(hHistMenu.Children);
        for i = startAt:endAt
            if isKey(obj.hClust.history, i)
                label = obj.hClust.history(i);
                uimenu(hHistMenu, 'Label', label, ...
                       'Callback', @(hO, hE) obj.restoreHistory(i), ...
                       'Checked', jrclust.utils.ifEq(i == hClust.editPos, 'on', 'off'), ...
                       'Enable', 'on');
            end
        end
    end
end

