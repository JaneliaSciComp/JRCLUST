%--------------------------------------------------------------------------
function mr = get_wav_(tmrWav_clu, viSite, iClu, viShift)
    if nargin <4, viShift = 0; end

    mrWav1 = (tmrWav_clu(:,viSite, iClu));
    % mrWav1 = zscore_(mrWav1);

    % induce waveform shift
    cmr = cell(1, numel(viShift));
    for iShift = 1:numel(viShift)
        mr1 = shift_mr_(mrWav1, viShift(iShift));
        mr1 = bsxfun(@minus, mr1, mean(mr1,1));
        cmr{iShift} = mr1(:);
    end
    % mr = zscore_(cell2mat_(cmr));
    mr = cell2mat_(cmr);
end %func
