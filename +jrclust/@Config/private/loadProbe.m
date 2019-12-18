function probe = loadProbe(probeFile)
    %LOADPROBE Load a probe file
    try % first try to load probe from mat file
        pstr = load(probeFile, '-mat');
    catch
        pstr = jrclust.utils.mToStruct(probeFile);
    end

    probe = struct();

    % exactly one of 'channels' or 'siteMap'
    if isfield(pstr, 'channels')
        if isfield(pstr, 'siteMap') && ~jrclust.utils.isEqual(pstr.channels(:), pstr.siteMap(:))
            error('both channels and siteMap specified, but they are not equal');
        end

        probe.siteMap = pstr.channels;
    elseif isfield(pstr, 'siteMap')
        probe.siteMap = pstr.siteMap;
    else
        error('missing channels/siteMap field');
    end

    % exactly one of 'geometry' or 'siteLoc'
    if isfield(pstr, 'geometry')
        if isfield(pstr, 'siteLoc') && ~jrclust.utils.isEqual(pstr.geometry(:), pstr.siteLoc(:))
            error('both geometry and siteLoc specified, but they are not equal');
        end

        probe.siteLoc = pstr.geometry;
    elseif isfield(pstr, 'siteLoc')
        probe.siteLoc = pstr.siteLoc;
    else
        error('missing geometry/siteLoc field');
    end

    % exactly one of 'pad' or 'probePad'
    if isfield(pstr, 'pad')
        if isfield(pstr, 'probePad') && ~jrclust.utils.isEqual(pstr.pad(:), pstr.probePad(:))
            error('both pad and probePad specified, but they are not equal');
        end

        probe.probePad = pstr.pad;
    elseif isfield(pstr, 'probePad')
        probe.probePad = pstr.probePad;
    else
        error('missing pad/probePad field');
    end

    % infer single shank from absence of shankMap
    if ~(isfield(pstr, 'shankMap') || isfield(pstr, 'shank') || isfield(pstr, 'cviShank'))
        probe.shankMap = ones(size(probe.siteMap));
    else
        if isfield(pstr, 'cviShank') && iscell(pstr, 'cviShank')
            pstr.cviShank = cell2mat(arrayfun(@(i) i*ones(size(pstr.cviShank{i})), 1:numel(pstr.cviShank), ...
                                     'UniformOutput', 0));
        end

        if isfield(pstr, 'shankMap') && isfield(pstr, 'cviShank')
            if ~jrclust.utils.isEqual(pstr.shankMap(:), pstr.cviShank(:))
                error('shankMap and cviShank specified, but they are not equal');
            end

            pstr = rmfield(pstr, 'cviShank');
        end

        if isfield(pstr, 'shankMap') && isfield(pstr, 'shank')
            if ~jrclust.utils.isEqual(pstr.shankMap(:), pstr.shank(:))
                error('shankMap and shank specified, but they are not equal');
            end

            pstr = rmfield(pstr, 'shank');
        end

        if isfield(pstr, 'shank') && isfield(pstr, 'cviShank')
            if ~jrclust.utils.isEqual(pstr.shank(:), pstr.cviShank(:))
                error('shank and cviShank specified, but they are not equal');
            end

            pstr.shankMap = pstr.shank;
            pstr = rmfield(pstr, 'shank');
            pstr = rmfield(pstr, 'cviShank');
        elseif isfield(pstr, 'shank')
            pstr.shankMap = pstr.shank;
        elseif isfield(pstr, 'cviShank')
            pstr.shankMap = pstr.cviShank;
        end

        probe.shankMap = pstr.shankMap;
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

