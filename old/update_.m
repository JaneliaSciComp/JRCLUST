%--------------------------------------------------------------------------
function update_(vcFile)
    % update_(vcFile) %update specific files
    % update_() %update files in default.cfg\sync_list

    S_cfg = read_cfg_();
    vcSource = S_cfg.path_dropbox;
    sprintf('copyfile ''%s\\%s'' .\\ f;', vcSource, 'default.cfg'); % copy default
    pause(.5); S_cfg = read_cfg_(); % reload default
    if strcmpi(pwd(), S_cfg.path_alpha), disp('cannot overwrite alpha'); return; end
    vcBackup = S_cfg.path_backup;
    mkdir_(vcBackup);

    t1 = tic;
    fprintf('Copying files to a backup location ''%s''...\n', vcBackup);
    copyfile_(fullfile('.', '*'), vcBackup); %backup

    fprintf('Updating from %s...\n', vcSource);
    copyfile_(S_cfg.sync_list, '.', vcSource);

    % Compile CUDA code
    fCompile = isempty(vcFile);
    try
        if fCompile, compile_cuda_(); end
    catch
        fprintf(2, 'CUDA code compilation error.\n');
    end

    fprintf('Updated, took %0.1fs.', toc(t1));
    fprintf('\tPrevious files backed up to %s\n', vcBackup);
    edit change_log.txt
end %func
