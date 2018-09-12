%--------------------------------------------------------------------------
function vrVrms_site = mr2rms_(mnWav2, max_sample)
    % uses median to estimate RMS
    if nargin<2, max_sample = []; end
    if isempty(max_sample)
        vrVrms_site = median(abs(mnWav2));
    else
        vrVrms_site = median(abs(subsample_mr_(mnWav2, max_sample, 1)));
    end
    vrVrms_site = single(vrVrms_site) / 0.6745;
end
