%--------------------------------------------------------------------------
function export_imec_sync_(vcFiles_prm)
    % sync output for IMEC (export the last channel entries)
    try
        P = file2struct_(vcFiles_prm);
        vcFile_bin = P.vcFile;
        fid = fopen(vcFile_bin, 'r');
        vnSync = fread(fid, inf, '*uint16'); fclose(fid);
        vnSync = vnSync(P.nChans:P.nChans:end); %subsample sync channel
        assignWorkspace_(vnSync);
    catch
        fprintf(2, 'error exporting sync: %s\n', lasterr());
    end
end
