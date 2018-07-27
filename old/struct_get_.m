%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function cvr = struct_get_(S, varargin)
    % Obtain a member of struct
    cvr = cell(size(varargin));
    for i=1:numel(varargin)
        vcName = varargin{i};
        if isfield(S, vcName)
            cvr{i} = S.(vcName);
        end
    end %for
end %func
