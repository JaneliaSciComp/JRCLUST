%--------------------------------------------------------------------------
function vrDist_k = knn_sorted_(mrFet_srt, n_neigh, k_nearest)
    nSpk = size(mrFet_srt, 1);
    % if nargin<3, n_neigh=[]; end
    % if nargin<4, k_nearest=[]; end

    vrDist_k = zeros(nSpk, 1, 'single');
    mrFet_srt = (mrFet_srt);
    n1 = n_neigh*2+1;
    mrFet_srt1 = mrFet_srt(1:n1,:);
    iCirc = 1;
    for iSpk = 1:nSpk
        %     iSpk = viSpk_sub(iSpk1);
        if iSpk > n_neigh && iSpk <= nSpk-n_neigh
            mrFet_srt1(iCirc,:) = mrFet_srt(iSpk+n_neigh,:);
            iCirc=iCirc+1;
            if iCirc>n1, iCirc=1; end
        end
        vrDist1 = sort(sum(bsxfun(@minus, mrFet_srt1, mrFet_srt(iSpk,:)).^2, 2));
        vrDist_k(iSpk) = vrDist1(k_nearest);
    end

    vrDist_k = gather_(vrDist_k);
end %func
