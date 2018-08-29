%--------------------------------------------------------------------------
function S = subsetStructElements(S, fieldNames, indices, dim)
    if isempty(fieldNames)
        return;
    end
    
    if nargin < 4
        dim = 1;
    end

    % function test
    if ischar(fieldNames)
        fieldNames = {fieldNames};
    end

    for i = 1:numel(fieldNames)
        fieldName = fieldNames{i};
        if ~isfield(S, fieldName)
            continue;
        end

        try
            val = S.(fieldName);
            if isempty(val)
                continue;
            end

            nValDims = ndims(val);
            if nValDims == 2 && any(size(val) == 1) % find a column or row vectors
                nValDims = 1;
            end

            switch nValDims
                case 1
                    val = val(indices);

                case 2
                    if dim == 1
                        val = val(indices, :);
                    elseif dim == 2
                        val = val(:, indices);
                    else
                        error('subsetStructElements: invalid dim');
                    end

                case 3
                    if dim == 1
                        val = val(indices, :, :);
                    elseif dim == 2
                        val = val(:, indices, :);
                    elseif dim == 3
                        val = val(:, :, indices);
                    else
                        error('subsetStructElements: invalid dim');
                    end
                otherwise
                    error('subsetStructElements: invalid # of dimensions (1-3 supported)');
            end % switch

            S.(fieldName) = val;
        catch
            error('subsetStructElements: %s field error', fieldName);
        end
    end % for
end % function
