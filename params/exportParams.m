%--------------------------------------------------------------------------
% 09/14/17 JJJ: Export settings to a prm file
function exportParams(paramFile, paramFileOut, fShow)
    % export current P structure to a parameter file
    % exportParams(paramFile): read and write to the paramFile (append default.prm)
    % exportParams(paramFile, paramFileOut): read from paramFile and output to paramFileOut
    if nargin < 3
        fShow = 1;
    end

    if isempty(paramFileOut)
        paramFileOut = paramFile;
    end

    copyfile(jrcpath_(read_cfg_('default_prm')), paramFileOut, 'f');

    P = get0_('P');
    if isempty(P)
        P = file2struct_(paramFile);
    end

    updateParamFile(P, paramFileOut);
    vcMsg = sprintf('Full parameter settings are exported to %s', paramFileOut);
    fprintf('%s\n', vcMsg);
    if fShow
        msgbox_(vcMsg);
        edit(paramFileOut);
    end
end %func
