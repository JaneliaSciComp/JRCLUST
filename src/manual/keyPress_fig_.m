%--------------------------------------------------------------------------
function keyPress_fig_(hFig, csKey)
    % Simulate key press function
    vcTag = get(hFig, 'Tag');
    % S0 = get(0, 'UserData');
    figure(hFig);
    figure_wait_(1);
    event1.Key = '';
    if ischar(csKey), csKey = {csKey}; end
    nKeys = numel(csKey);
    keyPressFcn__ = get(hFig, 'KeyPressFcn');
    for i=1:nKeys
        try
            event1.Key = csKey{i};
            keyPressFcn__(hFig, event1);
            fprintf('\tFigure ''%s'': Key ''%s'' success.\n', vcTag, csKey{i});
        catch
            fprintf(2, '\tFigure ''%s'': Key ''%s'' failed.\n', vcTag, csKey{i});
            disperr_();
        end
        %     pause(.1);
    end
    % drawnow;
    figure_wait_(0);
end %func
