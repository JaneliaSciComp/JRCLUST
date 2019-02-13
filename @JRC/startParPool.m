function startParPool(obj)
    %STARTPARPOOL Start the parallel pool
    if obj.hCfg.useParfor
        try
            parpool('local');
        catch % parpool already running
        end
    end
end