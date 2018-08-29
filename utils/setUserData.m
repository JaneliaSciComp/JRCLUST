%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function S0 = setUserData(varargin)
    % Set(0, 'UserData')

    S0 = get(0, 'UserData');
    % set(0, 'UserData', []); %prevent memory copy operation
    for i = 1:nargin
        try
            S0.(inputname(i)) = varargin{i};
        catch
            disperr_();
        end
    end
    set(0, 'UserData', S0);
end % function
