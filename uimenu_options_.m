%--------------------------------------------------------------------------
function menu_parent = uimenu_options_(menu_parent, csLabels, hFunc, hFig);
    % create options branch in the uimenu
    if nargin<4, hFig = []; end
    for i=1:numel(csLabels)
        uimenu(menu_parent, 'Label', csLabels{i}, 'Callback', @(h,e)hFunc(hFig, csLabels{i}, menu_parent));
    end
end % function
