%--------------------------------------------------------------------------
function vr = tnWav2uV_(vn, P)
    % Integrate the waveform

    if nargin<2, P = get0_('P'); end
    vr = jrclust.utils.bit2uV(vn, P);
    % if P.nDiff_filt>0,
    vr = jrclust.utils.meanSubtract(vr);
    % end
end %func
