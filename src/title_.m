%--------------------------------------------------------------------------
function hTitle = title_(hAx, vc)
    % title_(vc)
    % title_(hAx, vc)

    if nargin==1, vc=hAx; hAx=[]; end
    % Set figure title

    if isempty(hAx), hAx = gca; end
    hTitle = get_(hAx, 'Title');
    if isempty(hTitle)
        hTitle = title(hAx, vc, 'Interpreter', 'none', 'FontWeight', 'normal');
    else
        set_(hTitle, 'String', vc, 'Interpreter', 'none', 'FontWeight', 'normal');
    end
end %func
