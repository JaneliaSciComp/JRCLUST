%--------------------------------------------------------------------------
function vr = tnWav2uV_(vn, P)
    % Integrate the waveform

    if nargin<2, P = get0_('P'); end
    vr = bit2uV_(vn, P);
    % if P.nDiff_filt>0,
    vr = meanSubt_(vr);
    % end
end %func
