%--------------------------------------------------------------------------
function S = struct_reorder_(S, indices, varargin)
    for i = 1:numel(varargin)
        try
            fieldName = varargin{i};
            if ~isfield(S, fieldName)
                continue;
            end % ignore if not present

            fieldValue = S.(fieldName);

            if isvector(fieldValue)
                fieldValue = fieldValue(indices);
            elseif ismatrix(fieldValue)
                fieldValue = fieldValue(indices, :);
            else
                fieldValue = fieldValue(indices, :, :);
            end

            S.(fieldName) = fieldValue;
        catch
        end
    end
end %func
