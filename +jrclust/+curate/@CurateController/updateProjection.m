function updateProjection(obj, proj)
    %UPDATEPROJECTION Update the feature projection
    try
        obj.hCfg.dispFeature = proj;
    catch
        obj.hCfg.dispFeature = 'vpp';
    end

    % set menu items checked or unchecked where appropriate
    hProjMenu = obj.hMenus('ProjMenu');
    for i = 1:numel(hProjMenu.Children)
        if strcmp(hProjMenu.Children(i).Label, obj.hCfg.dispFeature)
            set(hProjMenu.Children(i), 'Checked', 'on');
        else
            set(hProjMenu.Children(i), 'Checked', 'off');
        end
    end

    obj.updateFigProj(1);
    obj.updateFigTime(1);
end
