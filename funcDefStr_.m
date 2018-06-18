%--------------------------------------------------------------------------
function P = funcDefStr_(P, varargin)
    csNames = varargin(1:2:end);
    csValues = varargin(2:2:end);

    for iField=1:numel(csNames)
        if ~isfield(P, csNames{iField})
            P = setfield(P, csNames{iField}, csValues{iField});
        end
    end
end %func
