function hClust = importv3(filename)
    %IMPORTV3
    hClust = [];

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

    hCfg = jrclust.Config2(prmFile);

    % construct dRes
    dRes = struct();
    errMsg = {};
    if ~isfield(S0, 'viTime_spk')
        errMsg{end+1} = 'missing viTime_spk';
    end
    if ~isfield(S0, 'viSite_spk')
        errMsg{end+1} = 'missing viTime_spk';
    end
    if ~isempty(errMsg)
        error(strjoin(errMsg, ';'));
    end
    dRes.spikeTimes = S0.viTime_spk;
    dRes.spikeSites = S0.viSite_spk;

    if isfield(S0, 'vrAmp_spk')
        dRes.spikeAmps = S0.vrAmp_spk;
    end
    if isfield(S0, 'mrPos_spk')
        dRes.spikePositions = S0.mrPos_spk;
    end
    if isfield(S0, 'cviSpk_site')
        dRes.spikesBySite = S0.cviSpk_site;
    end
    if isfield(S0, 'cviSpk2_site')
        dRes.spikesBySite2 = S0.cviSpk2_site;
    end
    if isfield(S0, 'cviSpk3_site')
        dRes.spikesBySite3 = S0.cviSpk3_site;
    end
    if isfield(S0, 'viSite2_spk')
        dRes.spikeSites2 = S0.viSite2_spk;
    end
    if isfield(S0, 'vrThresh_site')
        dRes.siteThresh = S0.vrThresh_site;
    end

    % load spkraw, spkwav, spkfet
    try
        fid = fopen(jrclust.utils.subsExt(prmFile, '_spkraw.jrc'), 'r');
        dRes.spikesRaw = fread(fid, inf, '*int16');
        dRes.spikesRaw = reshape(dRes.spikesRaw, S0.dimm_raw);
        fclose(fid);
    catch ME
        warning(ME.identifier, 'spikesRaw not imported: %s', ME.message);
    end

    try
        fid = fopen(jrclust.utils.subsExt(prmFile, '_spkwav.jrc'), 'r');
        dRes.spikesFilt = fread(fid, inf, '*int16');
        dRes.spikesFilt = reshape(dRes.spikesFilt, S0.dimm_spk);
        fclose(fid);
    catch ME
        warning(ME.identifier, 'spikesFilt not imported: %s', ME.message);
    end

    try
        fid = fopen(jrclust.utils.subsExt(prmFile, '_spkfet.jrc'), 'r');
        dRes.spikeFeatures = fread(fid, inf, '*single');
        dRes.spikeFeatures = reshape(dRes.spikeFeatures, S0.dimm_fet);
        fclose(fid);
    catch ME
        warning(ME.identifier, 'spikeFeatures not imported: %s', ME.message);
    end

    % construct sRes
    if isfield(S0, 'S_clu')
        S_clu = S0.S_clu;

        sRes = struct();
        reqFields = {'viClu', 'rho', 'delta', 'nneigh'};
        if all(ismember(reqFields, fieldnames(S_clu)))
            sRes.spikeClusters = S_clu.viClu;
            sRes.spikeRho = S_clu.rho;
            sRes.spikeDelta = S_clu.delta;
            sRes.spikeNeigh = S_clu.nneigh;
            if isfield(S_clu, 'ordrho')
                sRes.ordRho = S_clu.ordrho;
            else
                [~, sRes.ordRho] = sort(sRes.spikeRho, 'descend');
            end

            if isfield(S_clu, 'cviSpk_clu')
                sRes.spikesByCluster = S_clu.cviSpk_clu;
            end
            if isfield(S_clu, 'csNote_clu')
                sRes.clusterNotes = S_clu.csNote_clu;
            end
            if isfield(S_clu, 'icl')
                sRes.clusterCenters = S_clu.icl;
            end
            if isfield(S_clu, 'mrWavCor')
                sRes.simScore = S_clu.mrWavCor;
            end
            if isfield(S_clu, 'tmrWav_spk_clu')
                sRes.meanWfGlobal = S_clu.tmrWav_spk_clu;
            end
            if isfield(S_clu, 'tmrWav_raw_clu')
                sRes.meanWfGlobalRaw = S_clu.tmrWav_raw_clu;
            end
            if isfield(S_clu, 'trWav_spk_clu')
                sRes.meanWfLocal = S_clu.trWav_spk_clu;
            end
            if isfield(S_clu, 'trWav_raw_clu')
                sRes.meanWfLocalRaw = S_clu.trWav_raw_clu;
            end
            if isfield(S_clu, 'tmrWav_raw_hi_clu')
                sRes.meanWfRawHigh = S_clu.tmrWav_raw_hi_clu;
            end
            if isfield(S_clu, 'tmrWav_raw_hi_clu')
                sRes.meanWfRawLow = S_clu.tmrWav_raw_hi_clu;
            end
            if isfield(S_clu, 'viSite_clu')
                sRes.clusterSites = S_clu.viSite_clu;
            end
            if isfield(S_clu, 'viClu_auto')
                sRes.initialClustering = S_clu.viClu_auto;
            end
            if isfield(S_clu, 'viSite_min_clu')
                sRes.unitPeakSites = S_clu.viSite_min_clu;
            end
            if isfield(S_clu, 'vnSite_clu')
                sRes.nSitesOverThresh = S_clu.vnSite_clu;
            end
            if isfield(S_clu, 'vnSpk_clu')
                sRes.unitCount = S_clu.vnSpk_clu;
            end
            if isfield(S_clu, 'vrDc2_site')
                sRes.rhoCuts = S_clu.vrDc2_site;
            end
            if isfield(S_clu, 'vrIsiRatio_clu')
                sRes.unitISIRatio = S_clu.vrIsiRatio_clu;
            end
            if isfield(S_clu, 'vrIsoDist_clu')
                sRes.unitIsoDist = S_clu.vrIsoDist_clu;
            end
            if isfield(S_clu, 'vrLRatio_clu')
                sRes.unitLRatio = S_clu.vrLRatio_clu;
            end
            if isfield(S_clu, 'vrPosX_clu') && isfield(S_clu, 'vrPosY_clu')
                sRes.clusterCentroids = [S_clu.vrPosX_clu(:), S_clu.vrPosY_clu(:)];
            end
            if isfield(S_clu, 'vrSnr_clu')
                sRes.unitSNR = S_clu.vrSnr_clu;
            end
            if isfield(S_clu, 'vrVmin_clu')
                sRes.unitPeaks = S_clu.vrVmin_clu;
            end
            if isfield(S_clu, 'vrVmin_clu')
                sRes.unitPeaksRaw = S_clu.vrVmin_clu;
            end
            if isfield(S_clu, 'vrVpp_clu')
                sRes.unitVpp = S_clu.vrVpp_clu;
            end
            if isfield(S_clu, 'vrVpp_uv_clu')
                sRes.unitVppRaw = S_clu.vrVpp_uv_clu;
            end
            if isfield(S_clu, 'vrVrms_site')
                sRes.siteRMS = S_clu.vrVrms_site;
            end
        end
    end

    if ~isempty(sRes)
        hClust = jrclust.models.clustering.DensityPeakClustering(sRes, dRes, hCfg);
    else
        hClust = dRes;
    end
end