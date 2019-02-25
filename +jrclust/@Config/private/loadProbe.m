function probe = loadProbe(probeFile)
    %LOADPROBE Load a probe file
    try % first try to load probe from mat file
        pstr = load(probeFile, '-mat');
    catch
        pstr = jrclust.utils.mToStruct(probeFile);
    end

    probe = struct();

    assert(isfield(pstr, 'channels'), 'missing `channels` field');
    probe.siteMap = pstr.channels;

    assert(isfield(pstr, 'geometry'), 'missing `geometry` field');
    probe.siteLoc = pstr.geometry;

    assert(isfield(pstr, 'pad'), 'missing `pad` field');
    probe.probePad = pstr.pad;

    if ~isfield(pstr, 'shank')
        pstr.shank = [];
    end
    if ~isfield(pstr, 'cviShank')
        pstr.cviShank = [];
    end

    if isempty(pstr.shank) && ~isempty(pstr.cviShank)
        pstr.shank = pstr.cviShank;
    end
    if isempty(pstr.shank)
        probe.shankMap = ones(size(pstr.channels));
    elseif iscell(pstr.shank)
        probe.shankMap = cell2mat(arrayfun(@(i) i*ones(size(pstr.shank{i})), 1:numel(pstr.shank), 'UniformOutput', 0));
    else
        probe.shankMap = pstr.shank;
    end

    % load optional fields
    if isfield(pstr, 'maxSite')
        probe.nSiteDir = pstr.maxSite;
    end
    if isfield(pstr, 'nChans')
        probe.nChans = pstr.nChans;
    end
    if isfield(pstr, 'nSites_ref')
        probe.nSitesExcl = pstr.nSites_ref;
    end
    if isfield(pstr, 'sRateHz')
        probe.sampleRate = pstr.sRateHz;
    end
    if isfield(pstr, 'uV_per_bit')
        probe.bitScaling = pstr.uV_per_bit;
    end
    if isfield(pstr, 'um_per_pix')
        probe.umPerPix = pstr.um_per_pix;
    end
    if isfield(pstr, 'vcDataType')
        probe.dataType = pstr.vcDataType;
    end
    if isfield(pstr, 'viSiteZero')
        probe.ignoreSites = pstr.viSiteZero;
    end
end

