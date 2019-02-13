function mouseClickFigSim(obj, xyPos, clickType)
    %MOUSECLICKFIGSIM Handle callbacks for mouse clicks in sim view
    if obj.isWorking
        return;
    end

    xyPos = max(round(xyPos), [1 1]);
    if strcmp(clickType, 'normal') % left click
        obj.updateSelect(xyPos);
    end
end