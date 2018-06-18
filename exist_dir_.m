%--------------------------------------------------------------------------
% 11/5/17 JJJ: Created
function flag = exist_dir_(vcDir)
    flag = exist(vcDir, 'dir') == 7;
end %func
