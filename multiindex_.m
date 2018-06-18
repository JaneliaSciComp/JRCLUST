%--------------------------------------------------------------------------
function varargout = multiindex_(viKeep, varargin)
    % index first dimension of variables by a given index

    if nargout ~= numel(varargin), error('Number of argin=argout'); end

    if islogical(viKeep), viKeep = find(viKeep); end
    for i=1:numel(varargin)
        var1 = varargin{i};
        if isvector(var1)
            varargout{i} = var1(viKeep);
        elseif ismatrix(var1)
            varargout{i} = var1(viKeep,:);
        else %multidimensional variable
            varargout{i} = var1(viKeep,:,:);
        end
    end
end %func
