function success = autoMerge(obj)
%AUTOMERGE 
success = 0;
fprintf( "autoMerge for TemplateClustering\n");
if obj.nEdits > obj.editPos % not at tip of edit history, back out
    warning('cannot branch from history; use revert() first');
    return;
end

obj.orderClusters('clusterSites');

for iRepeat = 1:obj.hCfg.nPassesMerge % single-pass vs dual-pass correction
    nMerged = obj.mergeBySim();

    if nMerged < 1
        break;
    end
end
fprintf( 'size of waveform similarity matrix: \n')
size(obj.waveformSim)
fprintf( 'size of templateSim matrix: \n');
size(obj.templateSim)

success = 1;
end

