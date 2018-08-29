%--------------------------------------------------------------------------
function deleteMany(varargin)
    % provide cell or multiple arguments
    for i = 1:nargin
        try
            arg = varargin{i};

            if numel(arg) == 1
                delete(varargin{i});
            elseif iscell(arg)
                for j = 1:numel(arg)
                    try
                        delete(arg{j});
                    catch
                    end
                end
            else
                for j = 1:numel(arg)
                    try
                        delete(arg(j));
                    catch
                    end
                end
            end
        catch
        end
    end
end % function
