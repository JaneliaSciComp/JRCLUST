function tracesOut = bit2uV(tracesIn, hCfg)
    %BIT2UV Transform collected trace values to uV
    if strcmp(hCfg.filterType, 'sgdiff')
        nrm = sum((1:hCfg.nDiffOrder).^2) * 2;
    elseif strcmp(hCfg.filterType, 'ndiff')
        nrm = 2;
    else
        nrm = 1;
    end

    tracesOut = single(tracesIn) * single(hCfg.bitScaling / nrm);
end
