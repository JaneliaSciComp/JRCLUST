%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function import_tsf_(vcFile_tsf)
    % import tsf format (test spike file)
    % create a .bin file (fTranspose = 0)
    fprintf('Converting to WHISPER format (.bin and .meta)\n\t'); t1=tic;
    [mnWav, Sfile] = importTSF_(vcFile_tsf);
    % nChans = size(mnWav,2);
    vcFile_bin = strrep(vcFile_tsf, '.tsf', '.bin');
    write_bin_(vcFile_bin, mnWav);
    fprintf('\n\ttook %0.1fs.\n', toc(t1));

    % write .meta file
    vcFile_meta = jrclust.utils.subsExt(vcFile_tsf, '.meta');
    fid = fopen(vcFile_meta, 'W');
    fprintf(fid, 'niMNGain=200\n');
    fprintf(fid, 'niSampRate=%d\n', Sfile.sRateHz); %intan hardware default. in Smeta.header
    fprintf(fid, 'niAiRangeMax=0.6554\n'); %intan hardware default. in Smeta.header
    fprintf(fid, 'niAiRangeMin=-0.6554\n'); %intan hardware default. in Smeta.header
    fprintf(fid, 'nSavedChans=%d\n', Sfile.nChans); %intan hardware default. in Smeta.header
    fprintf(fid, 'fileTimeSecs=%f\n', Sfile.n_vd_samples/Sfile.sRateHz); %intan hardware default. in Smeta.header
    fprintf(fid, 'fileSizeBytes=%d\n', round(Sfile.n_vd_samples*2*Sfile.nChans)); %intan hardware default. in Smeta.header
    fclose(fid);

    % write meta file and bin file
    assignWorkspace_(Sfile);
    assignWorkspace_(mnWav);
    fprintf('Exported to %s, %s\n', vcFile_bin, vcFile_meta);
end %func
