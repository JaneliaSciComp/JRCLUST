%--------------------------------------------------------------------------
% 9/29/17 JJJ: Displaying the version number of the program and what's used. #Tested
function [vcVer, vcDate, vcVer_used] = jrcVersion(vcFile_prm)
    if nargin<1, vcFile_prm = ''; end
    vcVer = 'v3.2.5';
    vcDate = '1/8/2018';
    vcVer_used = '';
    if nargout==0
        fprintf('%s (%s) installed\n', vcVer, vcDate);
    end
    try
        if isempty(vcFile_prm)
            P = get0_('P');
            if ~isempty(P)
                fprintf('\t%s used in %s\n', P.version, P.prmFile);
                vcVer_used = P.version;
            end
        elseif exist_file_(vcFile_prm)
            P = loadParams(vcFile_prm);
            fprintf('\t%s used in %s\n', P.version, vcFile_prm);
            vcVer_used = P.version;
        end
    catch
        ;
    end
end %func
