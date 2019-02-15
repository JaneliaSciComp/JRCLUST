function success = exportToCSV(obj, zeroIndex)
    %EXPORTTOCSV Export spike times, sites, and clusters to CSV
    if nargin < 2
        zeroIndex = 0;
    end

    spikeTimes_ = double(obj.spikeTimes) / obj.hCfg.sampleRate;
    spikeSites_ = double(obj.spikeSites) - double(zeroIndex);

    filename = jrclust.utils.subsExt(obj.hCfg.configFile, '.csv');

    % write header
    try
        fid = fopen(filename, 'w');
        fprintf(fid, 'spikeTimes,spikeClusters,spikeSites\n');
        fclose(fid);
    catch ME
        warning('Failed to export: %s', ME.message);
        success = 0;
        return;
    end

    % write values
    dlmwrite(filename, [spikeTimes_(:), double(obj.spikeClusters(:)), spikeSites_(:)], 'precision', 9, '-append');

    if obj.hCfg.verbose
        fprintf('Wrote to %s. Columns:\n', filename);
        fprintf('\tColumn 1: Spike time (s)\n');
        fprintf('\tColumn 2: Unit# (positive #: valid units, 0: noise cluster, negative #: deleted clusters)\n');
        fprintf('\tColumn 3: Site# (starts with 1)\n');
    end

    success = 1;
end