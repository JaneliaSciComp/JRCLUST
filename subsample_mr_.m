%--------------------------------------------------------------------------
function [mr, vi] = subsample_mr_(mr, nMax, dimm)
    %[mr, vi] = subsample_mr_(mr, nMax, dimm)
    % subsample the column
    if nargin<3, dimm = 2; end
    if isempty(nMax), return ;end

    n = size(mr,dimm);
    nSkip = max(floor(n / nMax), 1);
    vi = 1:nSkip:n;
    if nSkip==1, return; end
    vi = vi(1:nMax);

    switch dimm
        case 2
        mr = mr(:,vi);
        case 1
        mr = mr(vi,:);
    end

    if nargout>=2
        if n > nMax
            vi = 1:nSkip:n;
            vi = vi(1:nMax);
        else
            vi = 1:n;
        end
    end
end % function
