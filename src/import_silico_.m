%--------------------------------------------------------------------------
function import_silico_(vcFile_prm, fSort)
    % need _gt struct?
    % import silico ground truth
    % S_gt: must contain viTime, viSite, viClu
    % [usage]
    % import_silico_(vcFile_prm, 0): use imported sort result (default)
    % import_silico_(vcFile_prm, 1): use jrclust sort result
    if nargin<2, fSort = 0; end

    global tnWav_raw tnWav_spk trFet_spk
    % convert jrc1 format (_clu and _evt) to jrc format. no overwriting
    % receive spike location, time and cluster number. the rest should be taken care by jrc processing
    P = loadParam_(vcFile_prm); %makeParam_kilosort_
    if isempty(P), return; end

    S_gt = load(strrep(vcFile_prm, '.prm', '_gt.mat')); %must contain viTime, viSite, viClu
    % vnSpk = cellfun(@numel, S.a);
    % viClu = int32(cell2mat_(arrayfun(@(n)n*ones(1, vnSpk(n)), 1:numel(vnSpk), 'UniformOutput', 0)));
    % viTime = int32(cell2mat_(S.a) * 20); % Convert to sample # (saved in ms unit & sampling rate =20KHZ)
    if ~isfield(S_gt, 'viSite')
        S_gt.viSite = S_gt.viSite_clu(S_gt.viClu);
    end
    S0 = struct('viTime_spk', S_gt.viTime(:), 'viSite_spk', S_gt.viSite(:), 'P', P, 'S_gt', S_gt);

    [tnWav_raw, tnWav_spk, trFet_spk, S0] = file2spk_(P, S0.viTime_spk, S0.viSite_spk);
    set(0, 'UserData', S0);

    % Save to file
    write_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), tnWav_raw);
    write_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), tnWav_spk);
    write_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), trFet_spk);

    % cluster and describe
    S0 = sort_(P);
    if ~fSort %use ground truth cluster
        S0.S_clu = S_clu_new_(S_gt.viClu, S0);
    end
    set(0, 'UserData', S0);

    % Save
    save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
    describe_(S0);
end %func
