function vcFile = default_prm_path_(failIfMissing)

    if nargin < 1
        failIfMissing = 1;
    end
    
    basedir = fileparts(jrcpath_());
    vcFile = fullfile(basedir, read_cfg_('default_prm'));
    
    if failIfMissing && exist(vcFile, 'file') ~= 2
        error('default param file %s is missing');
    end 
end
