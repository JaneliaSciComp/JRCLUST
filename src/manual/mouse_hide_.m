%--------------------------------------------------------------------------
function hidePlots = mouse_hide_(hFig, hObj_hide, S_fig)
    % hide during mouse pan to speed up
    if nargin<3, S_fig = get(hFig, 'UserData'); end
    % if nargin<3, S0 = get(0, 'UserData'); end
    if nargin == 0 %clear field
        %     try S_fig = rmfield(S_fig, 'vhFig_mouse'); catch; end
        try S_fig = rmfield(S_fig, 'hidePlots'); catch; end
    else
        if ~isfield(S_fig, 'vhFig_mouse') && ~isfield(S_fig, 'hidePlots')
            %         S_fig.vhFig_mouse = hFig;
            S_fig.hidePlots = {hObj_hide};
        else
            %         S_fig.vhFig_mouse(end+1) = hFig;
            S_fig.hidePlots{end+1} = hObj_hide;
        end
    end
    hidePlots = S_fig.hidePlots;
    if nargout==0, set(hFig, 'UserData', S_fig); end
end %func
