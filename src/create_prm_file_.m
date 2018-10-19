%--------------------------------------------------------------------------
% 7/31/17 JJJ: Handle wild-card in makeprm case
% 7/25/17 Create a parameter file from a meta file
% Probe file can be interpreted from SpikeGLX meta file.
function [P, vcPrompt] = create_prm_file_(vcFile_bin, vcFile_prb, vcFile_template, fAsk)
    if nargin<2, vcFile_prb = ''; end
    if nargin<3, vcFile_template = ''; end
    if nargin<4, fAsk = 1; end

    [P, vcPrompt] = deal([]);
    P0 = file2struct_(jrcpath_(read_cfg_('default_prm'))); %P = defaultParam();
    if ~isempty(vcFile_template)
        if exist(vcFile_template, 'file') == 2
            P0 = struct_merge_(P0, file2struct_(vcFile_template));
        else
            vcPrompt = sprintf('%s does not exist.\n', vcFile_bin);
            fprintf(2, '%s\n', vcPrompt);
            return;
        end
    end

    if any(vcFile_bin=='*') %wild card provided
        P.csFile_merge = vcFile_bin;
        vcFile_bin = strrep(vcFile_bin, '*', '');
    elseif isTextFile_(vcFile_bin)
        P.csFile_merge = vcFile_bin;
        %     vcFile_bin = subsFileExt_(vcFile_bin, '.bin');
    else
        if ~exist_file_(vcFile_bin)
            vcFile_bin_ = jrcpath_(vcFile_bin);
            if exist(vcFile_bin_, 'file') == 2, vcFile_bin = vcFile_bin_; end
        end
        if exist_file_(vcFile_bin)
            P.vcFile = vcFile_bin;
            P.csFile_merge = {};
        else
            vcPrompt = sprintf('%s does not exist.\n', vcFile_bin);
            fprintf(2, '%s\n', vcPrompt);
            return;
        end
    end

    % Load meta file
    if isempty(P.csFile_merge)
        vcFile_meta = subsFileExt_(vcFile_bin, '.meta');
    else
        csFiles_bin = filter_files_(P.csFile_merge);
        if isempty(csFiles_bin)
            vcFile_meta = '';
        else
            vcFile_meta = subsFileExt_(csFiles_bin{1}, '.meta');
        end
    end
    vcFile_meta = jrcpath_(vcFile_meta, 1);
    P_meta = read_meta_file_(vcFile_meta);
    if isempty(P_meta), P=[]; return; end

    % Get the probe file if missing
    if isempty(vcFile_prb)
        if isfield(P_meta.Smeta, 'imProbeOpt')
            if P_meta.Smeta.imProbeOpt > 0
                vcFile_prb = sprintf('imec3_opt%d.prb', round(P_meta.Smeta.imProbeOpt));
            end
        end
    end
    if isempty(vcFile_prb) % ask user
        vcFile_prb = inputdlg({'Probe file'}, 'Please specify a probe file', 1, {''});
        if isempty(vcFile_prb), P=[]; return; end
        vcFile_prb = vcFile_prb{1};
        if isempty(vcFile_prb)
            P = [];
            fprintf(2, 'You must specify a probe file (.prb)\n');
            return;
        end
    end

    % Assign prm file name
    [~,vcPostfix,~] = fileparts(vcFile_prb);
    P.vcFile_prm = subsFileExt_(vcFile_bin, ['_', vcPostfix, '.prm']);
    P.probe_file = vcFile_prb;
    try
        S_prb = file2struct_(find_prb_(vcFile_prb));
        %     P = struct_merge_(P, S_prb);
        if isfield(S_prb, 'maxSite'), P.maxSite = S_prb.maxSite; end
        if isfield(S_prb, 'nSites_ref'), P.nSites_ref = S_prb.nSites_ref; end
    catch
        disperr_(sprintf('Error loading the probe file: %s\n', vcFile_prb));
    end

    if exist(P.vcFile_prm, 'file') && fAsk
        vcAns = questdlg_('File already exists. Overwrite prm file?', 'Warning', 'Yes', 'No', 'No');
        if ~strcmpi(vcAns, 'Yes')
            P = [];
            vcPrompt = 'Cancelled by user.';
            return;
        end
    end

    % Load prb file
    if isfield(P, 'template_file')
        P = struct_merge_(file2struct_(P.template_file), P);
    end
    P = struct_merge_(P0, P);
    P = struct_merge_(P, P_meta);
    P = struct_merge_(P, file_info_(vcFile_bin));
    P.duration_file = P.nBytes_file / bytesPerSample_(P.vcDataType) / P.nChans / P.sRateHz; %assuming int16
    P.version = jrc_version_();
    try
        copyfile(jrcpath_(read_cfg_('default_prm')), P.vcFile_prm, 'f');
    catch
        fprintf(2, 'Invalid path: %s\n', P.vcFile_prm);
        return;
    end

    % Write to prm file
    edit_prm_file_(P, P.vcFile_prm);
    vcPrompt = sprintf('Created a new parameter file\n\t%s', P.vcFile_prm);
    disp(vcPrompt);
    if fAsk, edit(P.vcFile_prm); end % Show settings file
end %func
