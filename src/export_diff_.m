%--------------------------------------------------------------------------
function export_diff_(P)
    % export to _diff.bin, _diff.prb, _diff.prm files
    error('not implemented yet');
    if ~P.fTranspose_bin
        mnWav1 = reshape(load_bin_(P.vcFile, P.vcDataType), [], P.nChans);
        mnWav1 = mnWav1(:,P.viSite2Chan);
    else
        mnWav1 = reshape(load_bin_(P.vcFile, P.vcDataType), P.nChans, []);
        mnWav1 = mnWav1(P.viSite2Chan, :)';
    end

    % mnWav1: nT x nSites

    % fields to update, copy and save
    P1 = P;
    P1.vcFile = strrep(P.vcFile, '.bin', '_diff.bin');
    P1.vcFile_prm = strrep(P.vcFile_prm, '.prm', '_diff.prm');
    P1.probe_file = strrep(P.vcFile_prm, '.prm', '_diff.prb');
    P1.fTranspose_bin = 0;
    P1.vcCommonRef = 'none';
    P1.fDetectBipolar = 1;
    P1.nSites_ref = 0;

    % differentiate channels and write to bin file (two column type)
    nSites = numel(P.viSite2Chan);
    viChan_HP = 1:2:nSites;
    viChan_HN = 2:2:nSites;
    viChan_VP = 1:(nSites-4);
    viChan_VN = viChan_VP + 4;
    viChan_P = [viChan_HP(1), toRow_([viChan_VP(1:2:end); viChan_HP(2:end-1); viChan_VP(2:2:end)]), viChan_HP(end)];
    viChan_N = [viChan_HN(1), toRow_([viChan_VN(1:2:end); viChan_HN(2:end-1); viChan_VN(2:2:end)]), viChan_HN(end)];
    mnWav2 = mnWav1(:,viChan_P) - mnWav1(:,viChan_N);
    P1.nChans = size(mnWav2, 2);

    % Output files
    copyfile(P.vcFile_prm, P1.vcFile_prm, 'f');
    edit_prm_file_(P1, P1.vcFile_prm);
    write_bin_(P1.vcFile, mnWav2);
    % write to probe file
    % P1.probe_file
    % mnWav2 = load_bin_(strrep(P.vcFile, '.bin', '_diff.bin'), P.vcDataType); %read back test


end %func
