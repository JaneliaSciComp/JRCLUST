%--------------------------------------------------------------------------
% 2017/12/5 JJJ: distance-based neighboring unit selection
function mrWavCor = S_clu_wavcor_(S_clu, P, viClu_update)
    % version selector
    fUse_old = 1;

    spkLim_raw_factor = getOr(P, 'spkLim_raw_factor', 2);
    spkLim_factor_merge = getOr(P, 'spkLim_factor_merge', 1);
    if nargin<3, viClu_update = []; end

    if spkLim_factor_merge < spkLim_raw_factor && ~fUse_old
        mrWavCor = S_clu_wavcor_2_(S_clu, P, viClu_update); % newer format v3.1.8
    else
        mrWavCor = S_clu_wavcor_1_(S_clu, P, viClu_update); % older format
    end
end % function
