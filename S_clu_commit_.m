%--------------------------------------------------------------------------
function [S_clu, S0] = S_clu_commit_(S_clu, vcMsg)
    if nargin<2, vcMsg = ''; end
    if ~S_clu_valid_(S_clu)
        fprintf(2, '%s: Cluster number is inconsistent.', vcMsg);
        S0 = get(0, 'UserData');
        S_clu = get_(S0, 'S_clu');
    else
        S0 = setUserData(S_clu);
    end
end %func
