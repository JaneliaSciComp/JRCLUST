%--------------------------------------------------------------------------
% 8/9/17 JJJ: Generalized to any figure objects
function S0 = keyPressFcn_cell_(hObject, csKey, S0)
    % Simulate key press function

    if nargin<3, S0 = get(0, 'UserData'); end
    % figure_wait_(1);
    event1.Key = '';
    if ischar(csKey), csKey = {csKey}; end
    nKeys = numel(csKey);
    keyPressFcn_ = get(hObject, 'KeyPressFcn');
    for i=1:nKeys
        event1.Key = csKey{i};
        S0 = keyPressFcn_(hObject, event1, S0);
    end
    % drawnow;
    % figure_wait_(0);
    if nargout==0, set(0, 'UserData', S0); end
end %func
