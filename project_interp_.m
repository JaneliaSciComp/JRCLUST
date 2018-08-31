%--------------------------------------------------------------------------
function [mrFet1, mrFet2, mrFet3] = project_interp_(trWav2_spk, mrPv, P);
    [mrFet1, mrFet2, mrFet3] = deal([]);
    dimm1 = size(trWav2_spk);
    if ismatrix(trWav2_spk), dimm1(end+1) = 1; end
    mrWav_spk1 = reshape(trWav2_spk, dimm1(1), []);
    mrPv = gather_(mrPv);
    mrFet1 = reshape(mrPv(:,1)' * mrWav_spk1, dimm1(2:3))';
    if P.nPcPerChan >= 2
        mrFet2 = reshape(mrPv(:,2)' * mrWav_spk1, dimm1(2:3))';
    end
    if P.nPcPerChan >= 3
        mrFet3 = reshape(mrPv(:,3)' * mrWav_spk1, dimm1(2:3))';
    end

    % find optimal delay by interpolating 2x
    if ~getOr(P, 'fInterp_fet', 0), return; end
    vr1 = mrPv(:,1);
    vi0 = (1:numel(vr1))';
    viShift = [0, -1,-.5,.5,1]; %[0, -.5, .5]
    % viShift = [0, -1:.25:-.25,.25:.25:1];
    mrPv1 = zeros(numel(vr1), numel(viShift), 'like', mrPv);
    mrPv1(:,1) = vr1;
    for iShift = 2:numel(viShift)
        mrPv1(:,iShift) = zscore(interp1(vi0, vr1, vi0+viShift(iShift), 'pchip', 'extrap'));
    end
    mrPv1 = gpuArray_(mrPv1, isGpu_(trWav2_spk));
    [~, viMax_spk] = max(abs(mrPv1' * trWav2_spk(:,:,1)));
    for iShift=2:numel(viShift)
        viSpk2 = find(viMax_spk == iShift);
        if isempty(viSpk2), continue; end
        mrWav_spk2 = reshape(trWav2_spk(:,viSpk2,:), dimm1(1), []);
        mrFet1(:,viSpk2) = reshape(mrPv1(:,iShift)' * mrWav_spk2, [], dimm1(3))';
    end %for
end % function
