%--------------------------------------------------------------------------
function [S_clu, vlKeep_clu] = S_clu_remove_empty_(S_clu)

    vlKeep_clu = S_clu.vnSpk_clu>0;
    viClu_removed = find(~vlKeep_clu);
    if isempty(viClu_removed), return; end

    % waveform
    S_clu = S_clu_select_(S_clu, vlKeep_clu);
    if min(S_clu.viClu) < 1
        S_clu.viClu(S_clu.viClu<1) = 0;
        [~,~,S_clu.viClu] = unique(S_clu.viClu+1);
        S_clu.viClu = S_clu.viClu-1;
    else
        [~,~,S_clu.viClu] = unique(S_clu.viClu);
    end
    S_clu.viClu = int32(S_clu.viClu);
    S_clu.nClu = double(max(S_clu.viClu));
end %func
