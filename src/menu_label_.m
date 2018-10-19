%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function hMenu = menu_label_(vcTag, vcLabel)
    hMenu = [];
    try
        if ischar(vcTag)
            hMenu = findobj('Tag', vcTag, 'Type', 'uimenu');
        else
            hMenu = vcTag;
        end
        set(hMenu, 'Label', vcLabel); %figure property
    catch
        disperr_();
    end
end %func
