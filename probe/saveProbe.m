%--------------------------------------------------------------------------
% probes can now be saved to MAT files
function P = saveProbe(probeFilename, P)
    % append probe file to P
    if nargin < 2
        P = get0_('P');
    end

    % Find the probe file
    d = dir(probeFilename);
    if ~isempty(d) && d.bytes > 0
        answer = userDialog(sprintf('File ''%s'' exists and is not empty. Overwrite?', ...
            probeFilename));
        if ~strcmpi(answer, 'yes')
            return;
        end
    end

    [dirname, filename, ext] = fileparts(probeFilename);
    if isempty(dirname)
        probeFilename = fullfile(pwd, [filename, ext]);
    end

    P.probeFile = probeFilename;

    S_prb = struct('channels', P.chanMap, ...
                   'geometry', P.mrSiteXY, ...
                   'pad', P.vrSiteHW, ...
                   'shank', P.viShank_site);

    struct_save_(S_prb, probeFilename);
end %func
