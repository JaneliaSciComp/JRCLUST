%--------------------------------------------------------------------------
function S = struct_reorder_(S, viKeep, varargin)
    for i=1:numel(varargin)
        try
            vcVar = varargin{i};
            if ~isfield(S, vcVar), continue; end %ignore if not
            vr1 = S.(vcVar);
            if isvector(vr1)
                vr1 = vr1(viKeep);
            elseif ismatrix(vr1)
                vr1 = vr1(viKeep, :);
            else
                vr1 = vr1(viKeep, :, :);
            end
            S.(vcVar) = vr1;
        catch
            ;
        end
    end
end %func
