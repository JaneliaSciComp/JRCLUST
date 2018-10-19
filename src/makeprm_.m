%--------------------------------------------------------------------------
% 9/14/17 JJJ: supports custom template.prm file
% 8/2/17 JJJ: Testing and documentation
function [vcFile_prm, vcPrompt] = makeprm_(vcFile_bin, vcFile_prb, fAsk, vcFile_template)
    % Make a paramter file
    % vcFile_prm = makeprm_(vcFile_bin, vcFile_prb, fAsk, vcFile_template)
    % vcFile_prm = makeprm_(vcFile_bin, vcFile_template, fAsk)
    global fDebug_ui

    if nargin < 3
        fAsk = 1;
    end
    if nargin < 4
        vcFile_template = '';
    end

    if fDebug_ui == 1
        fAsk = 0;
    end
    [vcFile_prm, vcPrompt] = deal([]);

    if ~exist_file_(vcFile_bin)
        error('Data file %s does not exist\n', vcFile_bin);
    end

    set(0, 'UserData', []); % clear memory

    % given a parameter file instead of a probe
    if matchFileExt_(vcFile_prb, '.prm') && isempty(vcFile_template)
        vcFile_template = vcFile_prb;

        if ~exist_file_(vcFile_template)
            error('Template %s does not exist\n', vcFile_template);
        end

        S0 = file2struct_(vcFile_template);
        vcFile_prb = S0.probe_file;
    else % search for the probe
        if ~exist_file_(find_prb_(vcFile_prb))
            error('Probe %s does not exist\n', vcFile_prb);
        end
    end

    % data file, template, and probe found; make the parameter file
    [P, vcPrompt] = create_prm_file_(vcFile_bin, vcFile_prb, vcFile_template, fAsk);

    if isempty(P)
        [vcFile_prm, vcPrompt] = deal('');
        return;
    end

    set0_(P);
    vcFile_prm = P.vcFile_prm;
end
