%--------------------------------------------------------------------------
function export_lfp_(P)
    % export LFP waveform to workspace (ordered by the site numbers)

    P.vcFile_lfp = strrep(P.vcFile_prm, '.prm', '.lfp.jrc');
    if ~exist_file_(P.vcFile_lfp)
        import_lfp_(P)
    end
    mnLfp = load_bin_(P.vcFile_lfp, P.vcDataType);
    nSamples = floor(size(mnLfp,1) / P.nChans);
    mnLfp = reshape(mnLfp(1:P.nChans*nSamples), P.nChans, nSamples)';

    mnLfp = mnLfp(:, P.viSite2Chan);
    mrSiteXY = P.mrSiteXY;
    assignWorkspace_(mnLfp, mrSiteXY);
    fprintf('\tmnLfp has nSamples x nSites dimension, sites are ordered from the bottom to top, left to right\n');
    fprintf('\tmrSiteXY (site positions) has nSites x 2 dimension; col.1: x-coord, col.2: y-coord (um)\n');
end %func
