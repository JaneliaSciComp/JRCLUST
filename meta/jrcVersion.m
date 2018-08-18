%--------------------------------------------------------------------------
% 9/29/17 JJJ: Displaying the version number of the program and what's used. #Tested
function [versionStr, dateStr, versionUsed] = jrcVersion(paramFile)
    if nargin < 1
        paramFile = '';
    end

    versionStr = 'v3.2.5';
    dateStr = '1/8/2018';
    versionUsed = '';

    if nargout == 0
        fprintf('%s (%s) installed\n', versionStr, dateStr);
    end

    try
        if isempty(paramFile)
            P = get0_('P');
            if ~isempty(P)
                fprintf('\t%s used in %s\n', P.versionStr, P.paramFile);
                versionUsed = P.versionStr;
            end
        elseif fileExists(paramFile)
            P = loadParams(paramFile);
            fprintf('\t%s used in %s\n', P.versionStr, paramFile);
            versionUsed = P.versionStr;
        end
    catch
        ;
    end
end %func
