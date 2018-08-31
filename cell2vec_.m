%--------------------------------------------------------------------------
function vr = cell2vec_(cvr)
    % remove empty
    cvr = cvr(:);
    cvr = cvr(cellfun(@(x)~isempty(x), cvr));
    for i=1:numel(cvr)
        vr_ = cvr{i};
        cvr{i} = vr_(:);
    end
    vr = cell2mat(cvr);
end % function
