function md = info()
    %INFO Get JRCLUST repository metadata
    fid = fopen(fullfile(jrclust.utils.basedir(), 'json', 'jrc.json'), 'r', 'n', 'UTF-8');
    fstr = fread(fid, '*char')';
    fclose(fid);

    md = jsondecode(fstr);

    hash = gitHash();
    if ~isempty(hash)
        md.version.commitHash = hash;
    end
end

%% LOCAL FUNCTIONS
function h = gitHash()
    h = '';
    gitDir = fullfile(jrclust.utils.basedir(), '.git');
    if exist(gitDir, 'dir') == 7 && exist(fullfile(gitDir, 'HEAD'), 'file') == 2
        % get contents of HEAD
        fid = fopen(fullfile(gitDir, 'HEAD'), 'r', 'n', 'UTF-8');
        fstr = strip(fread(fid, '*char')');
        fclose(fid);

        refIdx = strfind(fstr, 'ref: '); % follow references
        if numel(refIdx) == 1 && refIdx == 1
            % remove 'ref: ' and use this path instead
            refPath = strsplit(strrep(fstr, 'ref: ', ''), '/');
            headPath = fullfile(gitDir, refPath{:});
            if exist(headPath, 'file') == 2
                fid = fopen(headPath, 'r', 'n', 'UTF-8');
                fstr = strip(fread(fid, '*char')');
                fclose(fid);
            end
        end

        if ~isempty(regexp(fstr, '^[0-9a-f]{5,40}$', 'once')) % probably a commit hash
            h = fstr;
        end
    end
end