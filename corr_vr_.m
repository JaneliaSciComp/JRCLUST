%--------------------------------------------------------------------------
function c = corr_vr_(vr1,vr2)
    % Vectorize matrix and compute correlation
    vr1 = vr1(:);
    vr2 = vr2(:);
    vr1 = vr1 - mean(vr1);
    vr2 = vr2 - mean(vr2);
    c = mean(vr1 .* vr2) / std(vr1,1) / std(vr2,1);
end %func
