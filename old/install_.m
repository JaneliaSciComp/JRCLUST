%--------------------------------------------------------------------------
% Install JRCLUST 1) create user.cfg 2) compile cuda
% 7/26/17 Code cleanup and test
function install_()
    % create user.cfg
    if exist('user.cfg', 'file') ~= 2
        fid = fopen('user.cfg', 'w');
        fprintf(fid, 'path_dropbox = ''C:\\Dropbox\\jrclust\\'';\n');
        fprintf(fid, 'path_backup = ''c:\\backup\\'';\n');
        fclose(fid);
        edit('user.cfg');
        msgbox_('Set path to ''path_dropbox'' and ''path_backup'' in user.cfg.');
    end

    compile_cuda_();
end %func
