%--------------------------------------------------------------------------
function [vlSuccess_menu, csLabel_menu] = menu_test_(hFig, csMenu_skip)
    vMenu0 = findobj('Type', 'uimenu', 'Parent', hFig);
    cvMenu = cell(size(vMenu0));
    for iMenu0 = 1:numel(vMenu0)
        cvMenu{iMenu0} = findobj('Type', 'uimenu', 'Parent', vMenu0(iMenu0))';
    end
    vMenu = [cvMenu{:}];
    cCallback_menu = get(vMenu, 'Callback');
    csLabel_menu = get(vMenu, 'Label');
    fprintf('\tTesting menu items\n');

    vlSuccess_menu = true(size(csLabel_menu));
    for iMenu = 1:numel(csLabel_menu)
        vcMenu = csLabel_menu{iMenu};
        if ismember(vcMenu, csMenu_skip), continue; end
        try
            hFunc = cCallback_menu{iMenu};
            if isempty(hFunc), continue; end
            %                     hFunc(hFigWav, []); %call function
            hFunc(vMenu(iMenu), []); %call function
            fprintf('\tMenu ''%s'' success.\n', vcMenu);
        catch
            fprintf(2, '\tMenu ''%s'' failed.\n', vcMenu);
            disperr_();
            vlSuccess_menu(iMenu) = 0;
        end
    end
end %func
