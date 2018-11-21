%--------------------------------------------------------------------------
function P = save_prb_(vcFile_prb, P)
    %SAVE_PRB_ save a probe to a MAT file

    % append probe file to P
    if nargin < 2
        P = get0_('P');
    end

    % Find the probe file
    d = dir(vcFile_prb);
    if ~isempty(d) && d.bytes > 0
        answer = questdlg_(sprintf('File ''%s'' exists and is not empty. Overwrite?', vcFile_prb));
        if ~strcmpi(answer, 'yes')
            return;
        end
    end

    [dirname, filename, ext] = fileparts(vcFile_prb);
    if isempty(dirname)
        vcFile_prb = fullfile(pwd, [filename, ext]);
    end

    P.probe_file = vcFile_prb;

    S_prb = struct('channels', P.viSite2Chan, ...
                   'geometry', P.mrSiteXY, ...
                   'pad', P.vrSiteHW, ...
                   'shank', P.viShank_site);

    jrclust.utils.saveStruct(S_prb, vcFile_prb);
end % function
