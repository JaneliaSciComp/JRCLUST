%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function hMenu = menu_checkbox_(vcTag, vcLabel, fUncheckRest)
    % Check the label selected and uncheck the rest
    if nargin<3, fUncheckRest = 1; end
    hMenu = [];
    try
        if ischar(vcTag)
            hMenu = findobj('Tag', vcTag, 'Type', 'uimenu');
        else
            hMenu = vcTag;
        end
        % Find children
        vhMenu_child = hMenu.Children;
        for iMenu = 1:numel(vhMenu_child)
            hMenu_child_ = vhMenu_child(iMenu);
            vcLabel_ = get(hMenu_child_, 'Label');
            if strcmpi(vcLabel_, vcLabel)
                set(hMenu_child_, 'Checked', 'on');
            elseif fUncheckRest
                set(hMenu_child_, 'Checked', 'off');
            end
        end %for
    catch
        disperr_();
    end
end %func
