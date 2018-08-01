%--------------------------------------------------------------------------
function reset_position_(hObject, event)
    % bottom to top left to right
    S0 = get(0, 'UserData');
    % P = S0.P;
    for iFig=1:numel(S0.figTags)
        hFig1 = getCachedFig(S0.figTags{iFig});
        if ishandle(hFig1)
            set(hFig1, 'OuterPosition', S0.figPositions{iFig});
            figure(hFig1); %bring it to foreground
        end
    end
end %func
