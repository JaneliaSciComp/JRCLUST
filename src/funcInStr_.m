%--------------------------------------------------------------------------
function P = funcInStr_( varargin )

    if isempty(varargin), P=struct(); return; end
    if isstruct(varargin{1}), P = varargin{1}; return; end

    csNames = varargin(1:2:end);
    csValues = varargin(2:2:end);
    P = struct();
    for iField=1:numel(csNames)
        if ~isfield(P, csNames{iField})
            %         v1 = csValues{iField};
            %         eval(sprintf('P.%s = v1;', csNames{iField}));
            P = setfield(P, csNames{iField}, csValues{iField});
        end
    end
end %func
