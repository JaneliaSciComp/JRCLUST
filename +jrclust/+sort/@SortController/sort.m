function res = sort(obj, dRes)
%SORT Cluster the spikes given in dRes
res = struct();
t0 = tic();

if ~isfield(dRes, 'spikeFeatures')
    obj.errMsg = 'cannot sort without features';
    obj.isError = 1;
    return;
end

nSpikes = numel(dRes.spikeTimes);

res.spikeRho = zeros(nSpikes, 1, 'single');
res.spikeDelta = zeros(nSpikes, 1, 'single');
res.spikeNeigh = zeros(nSpikes, 1, 'uint32');

% compute rho
res = jrclust.sort.computeRho(dRes, res, obj.hCfg);

% compute delta
res = jrclust.sort.computeDelta(dRes, res, obj.hCfg);

% assign clusters
[~, res.ordRho] = sort(res.spikeRho, 'descend');
res = jrclust.sort.assignClusters(dRes, res, obj.hCfg);
res.initialClustering = res.spikeClusters;

% create DensityPeakClusteringObject, compute derivate fields, automerge
hClust = jrclust.sort.DensityPeakClustering(obj.hCfg, res, dRes);
hClust.doRecompute();
hClust.autoMerge();

% don't recompute everything
if numel(hClust.recompute) > 0
    hClust.doRecompute();
end

res.hClust = hClust;

% if get_set_(P, 'fCorrect_overlap', 0) % correct waveforms and features after correcting clusters
%     S_clu = sort_overlap_(S0, S_clu, P);
% end

% summarize
res.sortTime = toc(t0);
res.sortedOn = now();
end % fun
