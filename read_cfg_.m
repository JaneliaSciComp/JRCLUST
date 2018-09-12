%--------------------------------------------------------------------------
function val = read_cfg_(vcName, fVerbose)
    % read configuration file that stores path to folder
    % load from default.cfg but override with user.cfg if it exists
    if nargin<2, fVerbose = 1; end

    S_cfg = file2struct_('default.cfg');
    % end
    if exist('user.cfg', 'file')
        S_cfg1 = file2struct_('user.cfg'); %override
        S_cfg = struct_merge_(S_cfg, S_cfg1, {'path_dropbox', 'path_backup', 'default_prm'});
        if fVerbose, fprintf('Configuration loaded from user.cfg.\n'); end
    else
        if fVerbose, fprintf('Configuration loaded from default.cfg.\n'); end
    end
    if nargin==0
        val = S_cfg;
    else
        try
            val = S_cfg.(vcName);
        catch
            disperr_('read_cfg_');
            switch lower(vcName)
                case 'default_prm'
                val = 'default.prm';
                otherwise
                val = [];
            end
        end
    end
end %func
