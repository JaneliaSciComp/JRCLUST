%--------------------------------------------------------------------------
% 9/17/17 JJJ: Created for SPAARC
function vcFile_prm = import_nsx_(vcFile_nsx, vcFile_prb, vcTemplate_prm)
    % Import neuroshare format
    % sample size is determined by the smallest file in the chan recording set
    if nargin<3, vcTemplate_prm = ''; end
    if matchFileExt_(vcFile_prb, '.prm')
        vcTemplate_prm = vcFile_prb;
        S_ = file2struct_(vcTemplate_prm);
        vcFile_prb = S_.probe_file;
    end

    % vcFile_nsx = 'E:\TimBruns\Ichabod Trial 14\exp_9_ichabod0014.ns5';
    if ~exist(vcFile_nsx, 'file'), error('File does not exist.') ;end
    [P, nSamples] = nsx_info_(vcFile_nsx);
    % [P, nSamples, vcFile_bin] = nsx2bin_(vcFile_nsx, 1);
    % P.fInverse_file = 1;
    % [mnWav, hFile, P] = load_nsx_(vcFile_nsx);
    P.probe_file = vcFile_prb;
    P.vcFile = vcFile_nsx;
    % mnWav = mnWav * -1; %inverse polarity
    [~, vcFile_prb_] = fileparts(vcFile_prb);
    vcFile_prm = subsFileExt_(P.vcFile, sprintf('_%s.prm', vcFile_prb_));
    if isempty(vcTemplate_prm)
        vcTemplate_prm = default_prm_path_();
    end
    assert_(exist_file_(vcTemplate_prm), sprintf('Template file does not exist: %s', vcTemplate_prm));

    % Write to a .prm file
    try
        S_prb = file2struct_(find_prb_(vcFile_prb));
        if isfield(S_prb, 'maxSite'), P.maxSite = S_prb.maxSite; end
        if isfield(S_prb, 'nSites_ref'), P.nSites_ref = S_prb.nSites_ref; end
    catch
        disperr_(sprintf('Error loading the probe file: %s\n', vcFile_prb));
    end
    P.duration_file = nSamples / P.sRateHz; %assuming int16
    P.version = jrclust.utils.version();
    P.vcFile_prm = vcFile_prm;
    % P.vcFile = vcFile_bin;
    copyfile(vcTemplate_prm, P.vcFile_prm, 'f');
    edit_prm_file_(P, P.vcFile_prm);
    vcPrompt = sprintf('Created a new parameter file\n\t%s', P.vcFile_prm);
    disp(vcPrompt);
    edit(P.vcFile_prm);
end
