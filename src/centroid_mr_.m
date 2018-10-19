%--------------------------------------------------------------------------
function [vrY, vrY2] = centroid_mr_(mrVpp, vrYe, mode1)
    % [vrX, vrY] = centroid_mr_(mrVpp, vrPos)
    % [mrXY] = centroid_mr_(mrVpp, vrPos)
    % mrVpp: nSites x nSpk

    if nargin<3, mode1 = 1; end
    if isrow(vrYe), vrYe = vrYe'; end
    % mrVpp_sq = abs(mrVpp);
    switch mode1
        case 1
        mrVpp = abs(mrVpp);
        case 2
        mrVpp = mrVpp.^2;
        case 3
        mrVpp = sqrt(abs(mrVpp));
        case 4
        mrVpp = mrVpp.^2;
        mrVpp = bsxfun(@minus, mrVpp, min(mrVpp));
    end

    vrVpp_sum = sum(mrVpp);
    % vrX = sum(bsxfun(@times, mrVpp_sq, mrSiteXY(:,1))) ./  vrVpp_sq_sum;
    vrY = sum(bsxfun(@times, mrVpp, vrYe)) ./  vrVpp_sum;
    if nargout>=2
        vrY2 = sum(bsxfun(@times, mrVpp, vrYe.^2)) ./  vrVpp_sum;
        vrY2 = sqrt(abs(vrY2 - vrY.^2));
    end
end %func
