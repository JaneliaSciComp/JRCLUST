%--------------------------------------------------------------------------
function csKeys = get_keyPress_(vcType)
    % return key press
    switch vcType
        case 'all'
        csKeys = [get_keyPress_('arrows'), get_keyPress_('misc'), get_keyPress_('alphanumeric')];
        case 'alphanumeric'
        csKeys = char([double('0'):double('9'), double('A'):double('Z'), double('a'):double('z')]);
        csKeys = num2cell(csKeys);
        case 'arrows'
        csKeys = {'uparrow', 'downarrow', 'leftarrow', 'rightarrow'};
        case 'misc'
        csKeys = {'home', 'end', 'space', 'esc'};
    end %switch
end %func
