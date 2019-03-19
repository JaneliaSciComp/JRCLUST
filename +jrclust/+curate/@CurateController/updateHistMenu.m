function updateHistMenu(obj)
%UPDATEHISTMENU Update the history menu
    if obj.hasMenu('HistMenu')
        hClust = obj.hClust;
        hHistMenu = obj.hMenus('HistMenu');

        endAt = hClust.nEdits + 1;
        startAt = max(1, endAt - 4); % show at most 5 (most recent) entries in hist menu

        delete(hHistMenu.Children);
        for i = startAt:endAt
            switch obj.hClust.history{i, 2}
                case 'delete'
                    label = sprintf('delete %s', jrclust.utils.field2str(obj.hClust.history{i, 3}));

                case 'merge'
                    label = sprintf('merge %d <- %d', obj.hClust.history{i, 3}, obj.hClust.history{i, 4});

                case 'split'
                    label = sprintf('split %d', obj.hClust.history{i, 3});

                case 'partition'
                    nSplits = numel(obj.hClust.history{i, 4}) - 1;
                    label = sprintf('partition %d (%d splits)', obj.hClust.history{i, 3}, nSplits);

                otherwise
                    label = obj.hClust.history{i, 2};
            end

            uimenu(hHistMenu, 'Label', label, ...
                   'Callback', @(hO, hE) obj.restoreHistory(i), ...
                   'Checked', jrclust.utils.ifEq(i - 1 == hClust.editPos, 'on', 'off'), ...
                   'Enable', 'on');
        end
    end
end

