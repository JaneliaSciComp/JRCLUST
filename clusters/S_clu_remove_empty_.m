%--------------------------------------------------------------------------
function [S_clu, vlKeep_clu] = S_clu_remove_empty_(S_clu, fSave_old)

    if nargin < 2, fSave_old = 0; end

    vlKeep_clu = S_clu.nSpikesPerCluster>0;
    if isempty(find(~vlKeep_clu, 1)), return; end

    % waveform
    S_clu = S_clu_select_(S_clu, vlKeep_clu);

    if fSave_old
        S_clu.spikeClustersOld = S_clu.spikeClusters; % save old 
    end

    if min(S_clu.spikeClusters) < 1
        S_clu.spikeClusters(S_clu.spikeClusters<1) = 0;
        [~,~,S_clu.spikeClusters] = unique(S_clu.spikeClusters+1);
        S_clu.spikeClusters = S_clu.spikeClusters-1;
    else
        [~,~,S_clu.spikeClusters] = unique(S_clu.spikeClusters);
    end
    S_clu.spikeClusters = int32(S_clu.spikeClusters);
    S_clu.nClusters = double(max(S_clu.spikeClusters));
end % function
