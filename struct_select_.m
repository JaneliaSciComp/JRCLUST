%--------------------------------------------------------------------------
function S = struct_select_(S, csNames, viKeep, iDimm)
    if isempty(csNames), return; end
    if nargin<4, iDimm = 1; end

    % function test
    % if nargin==0, S.a=rand(10,1);S.b=rand(10,3);S.c=rand(10,3,5); csNames={'a','b','c'}; viKeep=[1,2,3,5]; iDimm=1; end
    if ischar(csNames), csNames = {csNames}; end
    for i=1:numel(csNames)
        vcName_ = csNames{i};
        if ~isfield(S, vcName_), continue; end
        try
            val = S.(vcName_);
            if isempty(val), continue; end
            ndims_ = ndims(val);
            if ndims_==2 %find a column or row vectors
                if size(val,1)==1 || size(val,2)==1, ndims_=1; end %iscol or isrow
            end
            switch ndims_
                case 1, val = val(viKeep);
                case 2
                switch iDimm
                    case 1, val = val(viKeep,:);
                    case 2, val = val(:,viKeep);
                    otherwise
                    disperr_('struct_select_: invalid iDimm');
                end
                case 3
                switch iDimm
                    case 1, val = val(viKeep,:,:);
                    case 2, val = val(:,viKeep,:);
                    case 3, val = val(:,:,viKeep);
                    otherwise, disperr_('struct_select_: invalid iDimm');
                end
                otherwise, disperr_('struct_select_: invalid # of dimensions (1-3 supported)');
            end %switch
            S.(vcName_) = val;
        catch
            disperr_(sprintf('struct_select_: %s field error', vcName_));
        end
    end %for
end %func
