%--------------------------------------------------------------------------
function [S_clu, S0] = S_clu_new_(arg1, S0)
    % S_clu = S_clu_new_(S_clu, S0)
    % S_clu = S_clu_new_(arg1, S0)

    if nargin<2, S0 = get(0, 'UserData'); end
    % S_clu = S0.S_clu;
    % if ~isstruct(arg1), S_clu = struct(); end
    S_clu = get_(S0, 'S_clu'); %previous S_clu
    if isempty(S_clu), S_clu = struct(); end
    if ~isstruct(arg1)
        S_clu.viClu = arg1; %skip FigRD step for imported cluster
    else
        S_clu = struct_append_(S_clu, arg1);
    end

    S_clu.viClu = int32(S_clu.viClu);
    S_clu = S_clu_refresh_(S_clu);
    S_clu = S_clu_update_wav_(S_clu, S0.P);
    S_clu = S_clu_position_(S_clu);
    if ~isfield(S_clu, 'csNote_clu')
        S_clu.csNote_clu = cell(S_clu.nClusters, 1); %reset note
    end
    S_clu = S_clu_quality_(S_clu, S0.P);
    [S_clu, S0] = S_clu_commit_(S_clu, 'S_clu_new_');
end %func
