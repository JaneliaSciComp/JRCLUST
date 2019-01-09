function parent = menuOptions(parent, labels, hFun)
    %MENUOPTIONS Create options branch in the uimenu
    for i = 1:numel(labels)
        uimenu(parent, 'Label', labels{i}, 'Callback', @(hO, hE) hFun(labels{i}, parent));
    end
end