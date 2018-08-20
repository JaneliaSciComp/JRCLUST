%--------------------------------------------------------------------------
function export_diff_(P)
    % export to _diff.bin, _diff.prb, _diff.prm files
    error('not implemented yet');
    if ~P.fTranspose_bin
        mnWav1 = reshape(load_bin_(P.vcFile, P.dataType), [], P.nChans);
        mnWav1 = mnWav1(:,P.chanMap);
    else
        mnWav1 = reshape(load_bin_(P.vcFile, P.dataType), P.nChans, []);
        mnWav1 = mnWav1(P.chanMap, :)';
    end

    % mnWav1: nT x nSites

    % fields to update, copy and save
    P1 = P;
    P1.vcFile = strrep(P.vcFile, '.bin', '_diff.bin');
    P1.vcFile_prm = strrep(P.paramFile, '.prm', '_diff.prm');
    P1.probeFile = strrep(P.paramFile, '.prm', '_diff.prb');
    P1.fTranspose_bin = 0;
    P1.vcCommonRef = 'none';
    P1.fDetectBipolar = 1;
    P1.nSites_ref = 0;

    % differentiate channels and write to bin file (two column type)
    nSites = numel(P.chanMap);
    viChan_HP = 1:2:nSites;
    viChan_HN = 2:2:nSites;
    viChan_VP = 1:(nSites-4);
    viChan_VN = viChan_VP + 4;
    viChan_P = [viChan_HP(1), toRow_([viChan_VP(1:2:end); viChan_HP(2:end-1); viChan_VP(2:2:end)]), viChan_HP(end)];
    viChan_N = [viChan_HN(1), toRow_([viChan_VN(1:2:end); viChan_HN(2:end-1); viChan_VN(2:2:end)]), viChan_HN(end)];
    mnWav2 = mnWav1(:,viChan_P) - mnWav1(:,viChan_N);
    P1.nChans = size(mnWav2, 2);

    % Output files
    copyfile(P.paramFile, P1.vcFile_prm, 'f');
    updateParamFile(P1, P1.vcFile_prm);
    write_bin_(P1.vcFile, mnWav2);
    % write to probe file
    % P1.probeFile
    % mnWav2 = load_bin_(strrep(P.vcFile, '.bin', '_diff.bin'), P.dataType); %read back test


end %func
