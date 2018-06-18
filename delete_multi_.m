%--------------------------------------------------------------------------
function delete_multi_(varargin)
    % provide cell or multiple arguments
    for i=1:nargin
        try
            vr1 = varargin{i};
            if numel(vr1)==1
                delete(varargin{i});
            elseif iscell(vr1)
                for i1=1:numel(vr1)
                    try
                        delete(vr1{i1});
                    catch
                    end
                end
            else
                for i1=1:numel(vr1)
                    try
                        delete(vr1(i1));
                    catch
                    end
                end
            end
        catch
        end
    end
end %func
