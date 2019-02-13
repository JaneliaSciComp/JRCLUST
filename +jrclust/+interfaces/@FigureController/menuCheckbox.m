function hMenu = menuCheckbox(menuTag, label, doUncheckOthers)
    %MENUCHECKBOX Check the label selected, optionally uncheck the rest
    if nargin < 3
        doUncheckOthers = 1;
    end

    hMenu = [];
    try
        if ischar(menuTag)
            hMenu = findobj('Tag', menuTag, 'Type', 'uimenu');
        else
            hMenu = menuTag;
        end

        % Find children
        hMenuChildren = hMenu.Children;
        for iMenu = 1:numel(hMenuChildren)
            hMenuChild = hMenuChildren(iMenu);

            iLabel = get(hMenuChild, 'Label');
            if strcmpi(iLabel, label)
                set(hMenuChild, 'Checked', 'on');
            elseif doUncheckOthers
                set(hMenuChild, 'Checked', 'off');
            end
        end % for
    catch
    end
end