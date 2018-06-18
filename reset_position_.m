%--------------------------------------------------------------------------
function reset_position_(hObject, event)
    % bottom to top left to right
    S0 = get(0, 'UserData');
    % P = S0.P;
    for iFig=1:numel(S0.csFig)
        hFig1 = get_fig_cache_(S0.csFig{iFig});
        if ishandle(hFig1)
            set(hFig1, 'OuterPosition', S0.cvrFigPos0{iFig});
            figure(hFig1); %bring it to foreground
        end
    end
end %func
