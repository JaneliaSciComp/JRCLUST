%--------------------------------------------------------------------------
function commit_jrc3_(S_cfg, vcDir, fZipFile)
    global fDebug_ui

    if nargin<1, S_cfg = read_cfg_(); end
    if nargin<2, vcDir = []; end
    if nargin<3, fZipFile = 1; end

    fDebug_ui = 0;
    set0_(fDebug_ui);
    if isempty(vcDir)
        vcDir = [S_cfg.path_dropbox3, filesep()];
    end

    disp(['Commiting to ', vcDir]);
    for i=1:2
        disp(['Deleting files in ' vcDir]);
        delete_([vcDir, '*']);
        %     delete_([vcDir, '\kilosort\*']);
        pause(.5);
    end

    copyfile_(S_cfg.sync_list_ver3, vcDir); %destination root
    copyfile_(S_cfg.csFiles_cu3, vcDir); %destination root
    copyfile_(S_cfg.csFiles_ptx3, vcDir); %destination root
    % mkdir_([vcDir, 'kilosort']);
    % copyfile_('./kilosort/*', [vcDir, 'kilosort']); %destination root
    % copyfile_(S_cfg.csFiles_sample, vcDir);

    % Zip
    % delete_([vcDir, 'jrc3.zip']); %delete previously zipped file
    if fZipFile
        hMsg = msgbox_(sprintf('Archiving to %s', [vcDir, 'jrc3.zip']));
        t1 = tic;
        [csFiles_jrc3_full, csFiles_jrc3] = dir_([vcDir, '*'], 'jrc3.zip');
        zip([vcDir, 'jrc3.zip'], csFiles_jrc3, vcDir);
        fprintf('Zip file creation took %0.1f\n', toc(t1));
        close_(hMsg);
        msgbox_('Update the Dropbox link for www.jrclust.org');
    end
end %func
