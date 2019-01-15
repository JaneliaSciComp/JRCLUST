function probe = doLoadProbe(probeFile)
    %DOLOADPROBE Summary of this function goes here
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
        probe.shankMap = cell2mat(arrayfun(@(i) i*ones(size(pstr.shank{i})), 1:numel(shank), 'UniformOutput', 0));
    else
        probe.shankMap = pstr.shank;
    end
end

