%--------------------------------------------------------------------------
function [mrMin, mrMax] = getFet_spk_(viSpk1, viSites1, S0)
    % get feature for the spikes of interest

    if nargin<3, S0 = get(0, 'UserData'); end
    P = S0.P;

    switch lower(P.vcFet_show)
        case {'vmin', 'vpp'}
        tnWav_spk1 = tnWav2uV_(tnWav_spk_sites_(viSpk1, viSites1, S0), P);
        [mrMin, mrMax] = multifun_(@(x)abs(permute(x,[2,3,1])), min(tnWav_spk1), max(tnWav_spk1));
        case {'cov', 'spacetime'}
        [mrMin, mrMax] = calc_cov_spk_(viSpk1, viSites1);
        case 'pca'
        [mrMin, mrMax] = pca_pc_spk_(viSpk1, viSites1); %getall spikes whose center lies in certain range
        otherwise
        error('not implemented yet');
    end
end %func
