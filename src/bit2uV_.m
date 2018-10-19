%--------------------------------------------------------------------------
function vr_uV = bit2uV_(vn, P)
    % use only for filtered traces

    if nargin<2, P = get0_('P'); end
    % if isempty(get_(P, 'nDiff_filt')), P.nDiff_filt = 0; end
    switch lower(get_filter_(P))
        case 'sgdiff'
        norm = sum((1:P.nDiff_filt).^2) * 2;
        case 'ndiff'
        norm = 2;
        otherwise
        norm = 1;
    end
    % switch P.nDiff_filt
    %     case 0, norm = 1;
    %     case 1, norm = 2;
    %     case 2, norm = 10;
    %     case 3, norm = 28;
    %     otherwise, norm = 60;
    % end %switch
    vr_uV = single(vn) * single(P.uV_per_bit / norm);
end
