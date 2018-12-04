function [hCfg, dRes, sRes] = importv3(filename)
    %IMPORTV3
    [hCfg, dRes, sRes] = deal([]);

    dirname = fileparts(filename);
    if isempty(dirname)
        dirname = pwd();
    end

    try
        S0 = load(filename, '-mat');
    catch
        return;
    end

    % load config
    if ~isfield(S0, 'P') || ~isfield(S0.P, 'vcFile_prm')
        return;
    end
    prmFile = jrclust.utils.absPath(S0.P.vcFile_prm, dirname);
    if isempty(prmFile)
        return;
    end

    hCfg = jrclust.Config(prmFile);

    % construct dRes
    dRes = struct();
    dRes.spikeAmps = S0.vrAmp_spk;
    dRes.spikePositions = S0.mrPos_spk;
    dRes.spikesBySite = S0.cviSpk_site;
    dRes.spikesBySite2 = S0.cviSpk2_site;
    dRes.spikesBySite3 = S0.cviSpk3_site;
    dRes.spikeSites = S0.viSite_spk;
    dRes.spikeSites2 = S0.viSite2_spk;
    dRes.siteThresh = S0.vrThresh_site;
    dRes.spikeTimes = S0.viTime_spk;

    % load spkraw, spkwav, spkfet
    fid = fopen(strrep(prmFile, '.prm', '_spkraw.jrc'), 'r');
    dRes.spikesRaw = fread(fid, inf, '*int16');
    dRes.spikesRaw = reshape(dRes.spikesRaw, S0.dimm_raw);
    fclose(fid);

    fid = fopen(strrep(prmFile, '.prm', '_spkwav.jrc'), 'r');
    dRes.spikesFilt = fread(fid, inf, '*int16');
    dRes.spikesFilt = reshape(dRes.spikesFilt, S0.dimm_spk);
    fclose(fid);

    fid = fopen(strrep(prmFile, '.prm', '_spkfet.jrc'), 'r');
    dRes.spikeFeatures = fread(fid, inf, '*single');
    dRes.spikeFeatures = reshape(dRes.spikeFeatures, S0.dimm_fet);
    fclose(fid);

    % construct sRes
    sRes = struct();
    sRes.spikesByCluster = S0.S_clu.cviSpk_clu;
    sRes.clusterNotes = S0.S_clu.csNote_clu;
    sRes.spikeDelta = S0.S_clu.delta;
    sRes.clusterCenters = S0.S_clu.icl;
    sRes.simScore = S0.S_clu.mrWavCor;
    sRes.spikeNeigh = S0.S_clu.nneigh;
    sRes.ordRho = S0.S_clu.ordrho;
    sRes.spikeRho = S0.S_clu.rho;
    sRes.meanWfGlobal = S0.S_clu.tmrWav_spk_clu;
    sRes.meanWfGlobalRaw = S0.S_clu.tmrWav_raw_clu;
    sRes.meanWfLocal = S0.S_clu.trWav_spk_clu;
    sRes.meanWfLocalRaw = S0.S_clu.trWav_raw_clu;
    sRes.meanWfRawHigh = S0.S_clu.tmrWav_raw_hi_clu;
    sRes.meanWfRawLow = S0.S_clu.tmrWav_raw_hi_clu;
    sRes.spikeClusters = S0.S_clu.viClu;
    sRes.initialClustering = S0.S_clu.viClu_auto;
    sRes.clusterSites = S0.S_clu.viSite_clu;
    sRes.unitPeakSites = S0.S_clu.viSite_min_clu;
    sRes.nSitesOverThresh = S0.S_clu.vnSite_clu;
    sRes.clusterCounts = S0.S_clu.vnSpk_clu;
    sRes.unitISIRatio = S0.S_clu.vrIsiRatio_clu;
    sRes.unitIsoDist = S0.S_clu.vrIsoDist_clu;
    sRes.unitLRatio = S0.S_clu.vrLRatio_clu;
    sRes.clusterCentroids = [S0.S_clu.vrPosX_clu(:), S0.S_clu.vrPosY_clu(:)];
    sRes.unitSNR = S0.S_clu.vrSnr_clu;
    sRes.unitPeaks = S0.S_clu.vrVmin_clu;
    sRes.unitPeaksRaw = S0.S_clu.vrVmin_clu;
    sRes.unitVpp = S0.S_clu.vrVpp_clu;
    sRes.unitVppRaw = S0.S_clu.vrVpp_uv_clu;
    sRes.siteRMS = S0.S_clu.vrVrms_site;
end

% hCfg = jrclust.Config('F:\chen-20181126\SC009_112118_g0_t0.nidq_hh2_sc.prm');