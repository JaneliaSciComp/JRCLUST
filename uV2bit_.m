%--------------------------------------------------------------------------
function vr_uV = uV2bit_(vn, P)
    % use only for filtered traces

    if nargin<2, P = get0_('P'); end
    if ~isfield(P, 'nDiff_filt') || isempty(P.nDiff_filt), P.nDiff_filt = 0; end
    switch P.nDiff_filt
        case 0, norm = 1;
        case 1, norm = 2;
        case 2, norm = 10;
        case 3, norm = 28;
        otherwise, norm = 60;
    end %switch
    vr_uV = single(vn) / single(P.uV_per_bit / norm);
end
