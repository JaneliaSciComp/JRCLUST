%--------------------------------------------------------------------------
function [hFig, S_fig] = get_fig_(vcTag, hFig)
    % return figure handle based on the tag
    % cache figure handles
    % [usage]
    % get_fig_(vcTag)
    % get_fig_(vcTag, hFig) %build cache

    % multiple tags requested
    if iscell(vcTag), hFig = get_fig_all_(vcTag); return; end

    S_fig = [];
    try
        hFig = findobj('Tag', vcTag, 'Type', 'figure');
        if isempty(hFig)
            hFig = createFigure(vcTag);
            try %set position if exists
                S0 = get(0, 'UserData');
                iFig = find(strcmp(S0.csFig, vcTag), 1, 'first');
                if ~isempty(iFig)
                    set(hFig, 'OuterPosition', S0.cvrFigPos0{iFig});
                end
            catch
                %             disperr_();
            end
        else
            hFig=hFig(end); %get later one
        end
        if nargout>1, S_fig = get(hFig, 'UserData'); end
    catch
        hFig = [];
        disperr_();
    end
end %end
