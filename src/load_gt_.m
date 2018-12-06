%--------------------------------------------------------------------------
function S_gt = load_gt_(vcFile_gt, P)
    % S_gt contains viTime and viClu
    if nargin<2, P = get0_('P'); end
    if ~exist(vcFile_gt, 'file'), S_gt=[]; return; end
    S = load(vcFile_gt);
    if isfield(S, 'S_gt')
        S_gt = S.S_gt;
    elseif isfield(S, 'Sgt')
        S_gt = S.Sgt;
    elseif isfield(S, 'viClu') && isfield(S, 'viTime')
        S_gt = S;
    elseif isfield(S, 'viClu') && isfield(S, 'viSpk')
        S_gt.viTime = S.viSpk;
        S_gt.viClu = S.viClu;
    else
        % Convert Nick's format to JRCLUST fomat
        if isfield(S, 'gtTimes')
            S_gt.viTime = jrclust.utils.neCell2mat(S.gtTimes');
            S_gt.viClu = jrclust.utils.neCell2mat(arrayfun(@(i)ones(size(S.gtTimes{i}))*i, 1:numel(S.gtTimes), 'UniformOutput', 0)');
        else
            error('no field found.');
        end
        [S_gt.viTime, ix] = sort(S_gt.viTime, 'ascend');
        S_gt.viClu = S_gt.viClu(ix);
    end
    if ~isempty(get_(P, 'tlim_load'))
        nSamples = double(S_gt.viTime(end));
        nlim_load = min(max(round(P.tlim_load * P.sRateHz), 1), nSamples);
        viKeep = find(S_gt.viTime >= nlim_load(1) & S_gt.viTime <= nlim_load(2));
        [S_gt.viTime, S_gt.viClu] = multifun_(@(x)x(viKeep), S_gt.viTime, S_gt.viClu);
    end
    [viClu_unique, ~, viClu] = unique(S_gt.viClu);
    if max(S_gt.viClu) > numel(viClu_unique)
        S_gt.viClu = viClu;
    end
end %func
