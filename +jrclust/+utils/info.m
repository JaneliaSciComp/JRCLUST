function md = info()
    %INFO Get JRCLUST repository metadata
    fid = fopen(fullfile(jrclust.utils.basedir(), 'json', 'jrc.json'), 'r', 'n', 'UTF-8');
    fstr = fread(fid, '*char')';
    fclose(fid);

    md = jsondecode(fstr);
end

