%--------------------------------------------------------------------------
function disp_stats_(vr)
    % fprintf('n, mu/sd, (med,q25,q75), min-max:\t %d, %0.2f/%0.2f, (%0.2f, %0.2f, %0.2f), %0.2f-%0.2f\n', ...
    %     numel(vr), mean(vr), std(vr), quantile(vr, [.5, .25, .75]), min(vr), max(vr));
    fprintf('n, mu/sd, (10,25,*50,75,90%%), min-max:\t %d, %0.1f/%0.1f, (%0.1f, %0.1f, *%0.1f, %0.1f, %0.1f), %0.1f-%0.1f\n', ...
    numel(vr), mean(vr), std(vr), quantile(vr, [.1,.25,.5,.75,.9]), min(vr), max(vr));
end %func
