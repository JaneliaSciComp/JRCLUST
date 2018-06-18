%--------------------------------------------------------------------------
% 09/14/17 JJJ: Export settings to a prm file
function export_prm_(vcFile_prm, vcFile_out_prm, fShow)
    % export current P structure to a parameter file
    % export_prm_(vcFile_prm): read and write to the vcFile_prm (append default.prm)
    % export_prm_(vcFile_prm, vcFile_out_prm): read from vcFile_prm and output to vcFile_out_prm
    if nargin<3, fShow = 1; end

    if isempty(vcFile_out_prm), vcFile_out_prm = vcFile_prm; end
    copyfile(jrcpath_(read_cfg_('default_prm')), vcFile_out_prm, 'f');
    P = get0_('P');
    if isempty(P), P = file2struct_(vcFile_prm); end
    edit_prm_file_(P, vcFile_out_prm);
    vcMsg = sprintf('Full parameter settings are exported to %s', vcFile_out_prm);
    fprintf('%s\n', vcMsg);
    if fShow
        msgbox_(vcMsg);
        edit(vcFile_out_prm);
    end
end %func
