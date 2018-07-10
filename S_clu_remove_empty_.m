%--------------------------------------------------------------------------
function [S_clu, vlKeep_clu] = S_clu_remove_empty_(S_clu, fSave_old)

    if nargin < 2, fSave_old = 0; end

    vlKeep_clu = S_clu.vnSpk_clu>0;
    if isempty(find(~vlKeep_clu, 1)), return; end

    % waveform
    S_clu = S_clu_select_(S_clu, vlKeep_clu);
    
    if fSave_old
        S_clu.viClu_old = S_clu.viClu; % save old 
    end

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
