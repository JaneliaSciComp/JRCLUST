%--------------------------------------------------------------------------
% Use CPU
function [vl1, vii12] = find_overlap_spk_(viTime_clu1, viTime_clu12, nlimit, vrAmp_spk1, vrAmp_spk12)
    n1 = numel(viTime_clu1);
    n2 = numel(viTime_clu12);
    vii12 = zeros(size(viTime_clu1));
    i2prev = 1;
    vl1 = false(size(viTime_clu1));
    for iSpk1 = 1:numel(viTime_clu1)
        iTime1 = viTime_clu1(iSpk1);
        [vii_, i2prev] = findRange_(viTime_clu12, iTime1-nlimit, iTime1+nlimit, i2prev, n2);
        if isempty(vii_)
            continue;
        elseif numel(vii_) == 1
            if vrAmp_spk1(iSpk1) < vrAmp_spk12(vii_), continue; end
        else
            if any(vrAmp_spk1(iSpk1) < vrAmp_spk12(vii_)), continue; end
            [~, imin_] = min(abs(iTime1 - viTime_clu12(vii_)));
            vii_ = vii_(imin_);
        end
        vl1(iSpk1) = 1;
        vii12(iSpk1) = vii_;
    end %for
end %func
