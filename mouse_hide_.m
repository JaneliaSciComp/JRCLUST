%--------------------------------------------------------------------------
function cvhHide_mouse = mouse_hide_(hFig, hObj_hide, S_fig)
    % hide during mouse pan to speed up
    if nargin<3, S_fig = get(hFig, 'UserData'); end
    % if nargin<3, S0 = get(0, 'UserData'); end
    if nargin == 0 %clear field
        %     try S_fig = rmfield(S_fig, 'vhFig_mouse'); catch; end
        try S_fig = rmfield(S_fig, 'cvhHide_mouse'); catch; end
    else
        if ~isfield(S_fig, 'vhFig_mouse') && ~isfield(S_fig, 'cvhHide_mouse')
            %         S_fig.vhFig_mouse = hFig;
            S_fig.cvhHide_mouse = {hObj_hide};
        else
            %         S_fig.vhFig_mouse(end+1) = hFig;
            S_fig.cvhHide_mouse{end+1} = hObj_hide;
        end
    end
    cvhHide_mouse = S_fig.cvhHide_mouse;
    if nargout==0, set(hFig, 'UserData', S_fig); end
end % function
