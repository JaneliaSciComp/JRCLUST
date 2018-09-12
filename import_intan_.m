%--------------------------------------------------------------------------
% 9/17/17 JJJ: sample size is determined by the smallest file in the chan recording set
% 9/15/17 JJJ: Created
function vcFile_prm = import_intan_(vcFile_dat, vcFile_prb, vcArg3)

    % vcFile_dat = 'E:\ZGao\Exp_20170710_1302_run000\amp-A-000.dat';
    % vcFile_dat(end-7:end-5) = '000'; % replace to 000 index
    csFiles_dat = sort(dir_file_(vcFile_dat, 1)); % sort by the channel name
    [vcDir, vcFile] = fileparts(csFiles_dat{1});
    vcFile_bin = [vcDir, '.bin'];
    vcFile_prb = find_prb_(vcFile_prb);
    [~, vcFile_prb_] = fileparts(vcFile_prb);
    vcFile_prm = strrep(vcFile_bin, '.bin', sprintf('_%s.prm', vcFile_prb_));
    nChans = numel(csFiles_dat);
    nBytes_file = min(cellfun(@(vc)getBytes_(vc), csFiles_dat));
    P = struct('vcDataType', 'int16', 'probe_file', vcFile_prb, 'nChans', nChans, ...
    'uV_per_bit', .195, 'sRateHz', 30000, 'nBytes_file', nBytes_file);

    nSamples = P.nBytes_file / bytesPerSample_(P.vcDataType);

    % Read file and output
    mnWav  = zeros([nSamples, nChans], P.vcDataType);
    for iFile = 1:numel(csFiles_dat)
        try
            fid_ = fopen(csFiles_dat{iFile}, 'r');
            mnWav(iFile:nChans:end) = fread(fid_, inf, ['*', P.vcDataType]);
            fclose(fid_);
            fprintf('Loaded %s\n', csFiles_dat{iFile});
        catch
            fprintf(2, 'error %s\n', csFiles_dat{iFile});
        end
    end
    write_bin_(vcFile_bin, mnWav);
    clear mnWav;

    % Write to a .prm file
    try
        S_prb = file2struct_(vcFile_prb);
        if isfield(S_prb, 'maxSite'), P.maxSite = S_prb.maxSite; end
        if isfield(S_prb, 'nSites_ref'), P.nSites_ref = S_prb.nSites_ref; end
    catch
        disperr_(sprintf('Error loading the probe file: %s\n', vcFile_prb));
    end
    P.duration_file = nSamples / P.sRateHz; %assuming int16
    P.version = jrc_version_();
    P.vcFile_prm = vcFile_prm;
    P.vcFile = vcFile_bin;
    copyfile(jrcpath_(read_cfg_('default_prm')), P.vcFile_prm, 'f');
    edit_prm_file_(P, P.vcFile_prm);
    vcPrompt = sprintf('Created a new parameter file\n\t%s', P.vcFile_prm);
    disp(vcPrompt);
    edit(P.vcFile_prm);
end
