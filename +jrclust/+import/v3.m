function [hCfg, res] = v3(filename)
    %V3 Import an old-style session (_jrc.mat) to the new style
    [hCfg, res] = deal([]);

    filename_ = jrclust.utils.absPath(filename);
    if isempty(filename_)
        error('Could not find ''%s''', filename);
    end

    try
        S0 = load(filename_, '-mat');
    catch ME
        error('Failed to load ''%s'': %s', filename, ME.message);
    end

    % load config
    if ~isfield(S0, 'P')
        error('Could not find params struct (P)');
    end

    prmFile = jrclust.utils.absPath(S0.P.vcFile_prm, fileparts(filename_));
    if isempty(prmFile)
        % couldn't find vcFile_prm; set it manually
        dlgFieldNames = {'New config filename'};
        dlgFieldVals = {jrclust.utils.subsExt(strrep(filename_, '_jrc', ''), '.prm')};
        dlgAns = inputdlg(dlgFieldNames, 'Please specify a config filename', 1, dlgFieldVals, struct('Resize', 'on', 'Interpreter', 'tex'));
        if isempty(dlgAns) % abort
            return;
        end

        S0.P.vcFile_prm = dlgAns{1};
        % create the new config file
        [fid, errmsg] = fopen(S0.P.vcFile_prm, 'w');
        if fid == -1
            error('Could not open file ''%s'' for writing: %s', S0.P.vcFile_prm, errmsg);
        end
        fclose(fid);

        hCfg = jrclust.Config();
    else
        hCfg = jrclust.Config(prmFile);
    end

    % try to import params directly
    oldPrms = fieldnames(S0.P);

    % handle special cases in P
    if isfield(S0.P, 'spkLim_raw_ms') && isempty(S0.P.spkLim_raw_ms)
        if isfield(S0.P, 'spkLim_ms') && isfield(S0.P, 'spkLim_raw_factor')
            S0.P.spkLim_raw_ms = S0.P.spkLim_ms * S0.P.spkLim_raw_factor;
        end
    end

    for i = 1:numel(oldPrms)
        propname = oldPrms{i};

        if isfield(hCfg.oldParamSet, propname) && ~isempty(S0.P.(propname))
             % these will be converted automatically
             try
                hCfg.(propname) = S0.P.(propname);
             catch % use default
             end
        end
    end

    % construct dRes
    res = struct();
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
    res.spikeTimes = S0.viTime_spk;
    res.spikeSites = S0.viSite_spk;

    if isfield(S0, 'vrAmp_spk')
        res.spikeAmps = S0.vrAmp_spk;
    end
    if isfield(S0, 'mrPos_spk')
        res.spikePositions = S0.mrPos_spk;
    end
    if isfield(S0, 'cviSpk_site')
        res.spikesBySite = S0.cviSpk_site;
    end
    if isfield(S0, 'cviSpk2_site')
        res.spikesBySite2 = S0.cviSpk2_site;
    end
    if isfield(S0, 'cviSpk3_site')
        res.spikesBySite3 = S0.cviSpk3_site;
    end
    if isfield(S0, 'viSite2_spk')
        res.spikeSites2 = S0.viSite2_spk;
    end
    if isfield(S0, 'vrThresh_site')
        res.siteThresh = S0.vrThresh_site;
    end
    if isfield(S0, 'dimm_raw')
        res.rawShape = S0.dimm_raw;
    end
    if isfield(S0, 'dimm_spk')
        res.filtShape = S0.dimm_spk;
    end
    if isfield(S0, 'dimm_fet')
        res.featuresShape = S0.dimm_fet;
    end

    % rename spkraw, spkwav, spkfet
    spkraw = jrclust.utils.subsExt(strrep(filename_, '_jrc', ''), '_spkraw.jrc');
    if exist(spkraw, 'file') == 2
        renameFile(spkraw, jrclust.utils.subsExt(strrep(filename_, '_jrc', ''), '_raw.jrc'));
    else
        warning('Could not find file ''%s''', spkraw);
    end

    spkwav = jrclust.utils.subsExt(strrep(filename_, '_jrc', ''), '_spkwav.jrc');
    if exist(spkwav, 'file') == 2
        renameFile(spkwav, jrclust.utils.subsExt(strrep(filename_, '_jrc', ''), '_filt.jrc'));
    else
        warning('Could not find file ''%s''', spkwav);
    end

    spkfet = jrclust.utils.subsExt(strrep(filename_, '_jrc', ''), '_spkfet.jrc');
    if exist(spkfet, 'file') == 2
        renameFile(spkfet, jrclust.utils.subsExt(strrep(filename_, '_jrc', ''), '_features.jrc'));
    else
        warning('Could not find file ''%s''', spkfet);
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
            else
                sRes.spikesByCluster = arrayfun(@(iC) find(sRes.spikeClusters == iC), 1:max(sRes.spikeClusters), 'UniformOutput', 0);
            end
            if isfield(S_clu, 'csNote_clu')
                sRes.clusterNotes = S_clu.csNote_clu;
            end
            if isfield(S_clu, 'icl')
                clusterCenters = zeros(max(sRes.spikeClusters), 1);
                for iCluster = 1:numel(clusterCenters)
                    iSpikes = sRes.spikesByCluster{iCluster};
                    iCenters = intersect(S_clu.icl, iSpikes);
                    if isempty(iCenters)
                        [~, densest] = max(sRes.spikeRho(iSpikes));
                        clusterCenters(iCluster) = iSpikes(densest);
                    else
                        [~, densest] = max(sRes.spikeRho(iCenters));
                        clusterCenters(iCluster) = iCenters(densest);
                    end
                end
                sRes.clusterCenters = clusterCenters;
            end
            if isfield(S_clu, 'mrWavCor')
                sRes.waveformSim = S_clu.mrWavCor;
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
                sRes.rhoCutSite = S_clu.vrDc2_site;
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

        hClust = jrclust.sort.DensityPeakClustering(sRes, res, hCfg);
        msgs = hClust.inconsistentFields();
        assert(isempty(msgs), strjoin(msgs, ', '));

        % remove quality scores/initial clustering from sRes
        sRes = rmfield(sRes, {'clusterCentroids', 'clusterNotes', 'initialClustering', ...
                              'meanWfGlobal', 'meanWfGlobalRaw', 'meanWfLocal', 'meanWfLocalRaw', ...
                              'meanWfRawHigh', 'meanWfRawLow', 'nSitesOverThresh', 'waveformSim', ...
                              'siteRMS', 'unitISIRatio', 'unitIsoDist', 'unitLRatio', 'unitPeakSites', ...
                              'unitPeaks', 'unitPeaksRaw', 'unitSNR', 'unitVpp', 'unitVppRaw'});
        res = jrclust.utils.mergeStructs(res, sRes);
        res.hClust = hClust;
    end
end

%% LOCAL FUNCTIONS
function renameFile(oldfile, newfile)
    try
        movefile(oldfile, newfile);
    catch ME
        warning('failed to rename %s to %s: %s (try to move it manually?)', oldfile, newfile, ME.message);
    end
end