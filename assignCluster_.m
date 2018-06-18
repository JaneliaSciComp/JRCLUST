%--------------------------------------------------------------------------
function [cl, icl] = assignCluster_(cl, ordrho, nneigh, icl)
    ND = numel(ordrho);
    nClu = numel(icl);

    if isempty(cl)
        cl = zeros([ND, 1], 'int32');
        cl(icl) = 1:nClu;
    end

    if numel(icl) == 0 || numel(icl) == 1
        cl = ones([ND, 1], 'int32');
        icl = ordrho(1);
    else
        nneigh1 = nneigh(ordrho);
        for i=1:10
            vi = find(cl(ordrho)<=0);
            if isempty(vi), break; end
            vi=vi(:)';
            for ii = vi
                cl(ordrho(ii)) = cl(nneigh1(ii));
            end
            n1 = sum(cl<=0);
            if n1==0, break; end
            fprintf('i:%d, n0=%d, ', i, n1);
        end
        cl(cl<=0) = 1; %background
    end
end %func
