%--------------------------------------------------------------------------
% 9/19/17 JJJ: Created for SPARC
function mrRate_clu = clu_rate_(S_clu, viClu, nSamples)
    S0 = get(0, 'UserData');
    P = S0.P;
    if nargin<2, viClu = []; end
    if nargin<3, nSamples = []; end

    if isempty(viClu), viClu = 1:S_clu.nClusters; end
    sampleRateHz_rate = getOr(P, 'sampleRateHz_rate', 1000);
    filter_sec_rate = getOr(P, 'filter_sec_rate', 1);
    nFilt = round(sampleRateHz_rate * filter_sec_rate / 2);
    filter_shape_rate = lower(getOr(P, 'filter_shape_rate', 'triangle'));
    switch filter_shape_rate
        case 'triangle'
        vrFilt = ([1:nFilt, nFilt-1:-1:1]'/nFilt*2/filter_sec_rate);
        case 'rectangle'
        vrFilt = (ones(nFilt*2, 1) / filter_sec_rate);
        %     case 'interval3'
        %         1/2 1/3 1/6. interval history
    end %switch
    vrFilt = single(vrFilt);

    if isempty(nSamples), nSamples = round(P.duration_file * sampleRateHz_rate); end
    mrRate_clu = zeros([nSamples, numel(viClu)], 'single');
    for iClu1 = 1:numel(viClu)
        iClu = viClu(iClu1);
        viSpk_clu = S_clu.spikesByCluster{iClu};
        viTime_clu = S0.spikeTimes(viSpk_clu);
        viTime_ = round(double(viTime_clu) / P.sampleRateHz * sampleRateHz_rate);
        viTime_ = max(min(viTime_, nSamples), 1);
        mrRate_clu(viTime_, iClu1) = 1;
        mrRate_clu(:,iClu1) = conv(mrRate_clu(:,iClu1), vrFilt, 'same');
    end
end % function
