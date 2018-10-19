%--------------------------------------------------------------------------
function mi = rankorder_mr_(mr, val0)
    if nargin<2, val0 = 0; end % separate positive and negaitve number ranks

    if isrow(mr), mr=mr'; end

    dimm1 = size(mr); %=numel(vr);
    if numel(dimm1)==3
        %     mr = reshape(mr, [], dimm1(3));  % spike by spike order
        mr = mr(:); % global order
    end

    mi=zeros(size(mr));
    for iCol = 1:size(mr,2)
        vr1 = mr(:,iCol);
        vi_p = find(vr1>val0);
        if ~isempty(vi_p)
            [~, vi_p_srt] = sort(vr1(vi_p), 'ascend');
            mi(vi_p(vi_p_srt), iCol) = 1:numel(vi_p);
        end

        vi_n = find(vr1<val0);
        if ~isempty(vi_n)
            [~, vi_n_srt] = sort(vr1(vi_n), 'descend');
            mi(vi_n(vi_n_srt), iCol) = -(1:numel(vi_n));
        end
    end
    % [~,miSort] = sort(mr, vcOrder);
    % vr_ = (1:size(mr,1))';
    % for iCol = 1:size(mr,2)
    %     vi_ = miSort(:,iCol);
    %     mr(vi_, iCol)
    %     mi(,iCol) = vr_;
    % end

    if numel(dimm1)==3, mi = reshape(mi, dimm1); end
end %func
