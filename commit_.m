%--------------------------------------------------------------------------
% update jrclust_alpha to jrclust directory
function commit_(vcArg1)
    % commit_()
    %   Full validation and update
    % commit_('log')
    %   Update 'change_log.txt' only
    % commit_('skip')
    %   Skip unit test

    if nargin<1, vcArg1=''; end
    t1 = tic;
    S_cfg = read_cfg_();
    if ~strcmpi(pwd(), S_cfg.path_alpha), disp('must commit from alpha'); return; end

    if strcmpi(vcArg1, 'log')
        %     sprintf('copyfile change_log.txt ''%s'' f;', S_cfg.path_dropbox);
        copyfile_('change_log.txt', ...
        struct_get_(S_cfg, ...
        'path_dropbox', 'path_dropbox2', 'path_web', 'path_bitbucket', 'path_github'));
        disp('Commited change_log.txt');
        return;
    elseif ~strcmpi(vcArg1, 'skip')
        disp('Running unit tests before commit... ');
        nTests_fail = unit_test_(); %run full unit test
        if nTests_fail > 0
            fprintf(2, 'Commit aborted, %d unit tests failed\n.', nTests_fail);
            return;
        end
    else
        fprintf(2, 'Skipping unit test...\n');
    end

    % commit all
    disp(['Commiting to ', S_cfg.path_dropbox]);
    delete_files_(find_empty_files_());
    copyfile_(S_cfg.sync_list, S_cfg.path_dropbox); %destination root
    % copyfile_('./kilosort/*', [S_cfg.path_dropbox, filesep(), 'kilosort', filesep()]); %destination root
    copyfile_(S_cfg.csFiles_sample, S_cfg.path_dropbox);

    % commit jrc3 related files only
    try commit_jrc3_(S_cfg, [S_cfg.path_web, filesep()], 1); catch; disperr_(); end
    try commit_jrc3_(S_cfg, [S_cfg.path_bitbucket, filesep()], 0); catch; disperr_(); end
    try commit_jrc3_(S_cfg, [S_cfg.path_github, filesep()], 0); catch; disperr_(); end

    edit change_log.txt
    fprintf('Commited, took %0.1fs.\n', toc(t1));
end
