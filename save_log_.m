%--------------------------------------------------------------------------
function S0 = save_log_(vcCmd, S0)

    % save cluster info and save to file (append)
    % check for crash
    % todo: save differential increment from the beginning

    if nargin<2, S0 = get(0, 'UserData'); end
    [cS_log, P, S_clu, miClu_log] = get_(S0, 'cS_log', 'P', 'S_clu', 'miClu_log');

    if ~isempty(strfind(vcCmd, 'annotate'))
        S_log = cS_log{end};
        S_log.csNote_clu = S_clu.csNote_clu;
        cS_log{end} = S_log;
    else
        S_log = struct_('vcCmd', vcCmd, 'datenum', now(), 'csNote_clu', S_clu.csNote_clu);
        if isempty(cS_log) || ~iscell(cS_log)
            cS_log = {S_log};
        else
            cS_log{end+1} = S_log;
        end
    end

    % Keep P.MAX_LOG history
    if isempty(miClu_log)
        miClu_log = zeros([numel(S_clu.viClu), P.MAX_LOG], 'int16');
    end
    miClu_log(:, 2:end) = miClu_log(:, 1:end-1);
    miClu_log(:, 1) = int16(S_clu.viClu);
    %struct_save_(strrep(P.vcFile_prm, '.prm', '_log.mat'), 'cS_log', cS_log);
    S_log.viClu = int16(S_clu.viClu);
    struct_save_(S_log, strrep(P.vcFile_prm, '.prm', '_log.mat'), 0);
    S0.cS_log = cS_log;
    S0.miClu_log = miClu_log;

    ui_update_log_(cS_log, S0); % update revert to list

    if nargout<1, set(0, 'UserData', S0); end
end %func
