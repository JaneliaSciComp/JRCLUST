%--------------------------------------------------------------------------
function export_chan_(P, vcArg1)
    % export list of channels to a bin file, use fskip?
    if isempty(vcArg1)
        vcArg1 = inputdlg('Which channel(s) to export (separate by commas or space)', 'Channel', 1, {num2str(P.nChans)});
        if isempty(vcArg1), return; end
        vcArg1 = vcArg1{1};
    end
    viChan = str2num(vcArg1);
    if isnan(viChan), fprintf(2, 'Must provide a channel number\n'); return; end
    if any(viChan > P.nChans | viChan < 0)
        fprintf(2, 'Exceeding nChans (=%d).\n', P.nChans);
        return;
    end
    try
        vcChan_ = sprintf('%d-', viChan);
        vcFile_out = strrep(P.vcFile_prm, '.prm', sprintf('_ch%s.jrc', vcChan_(1:end-1)));
        mn = load_bin_chan_(P, viChan);
        write_bin_(vcFile_out, mn);
        if numel(viChan) == 1
            eval(sprintf('ch%d = mn;', viChan));
            eval(sprintf('assignWorkspace_(ch%d);', viChan));
        else
            assignWorkspace_(mn);
        end
    catch
        fprintf('Out of memory, exporting individual channels\n');
        for iChan1 = 1:numel(viChan)
            iChan = viChan(iChan1);
            vcFile_out = strrep(P.vcFile_prm, '.prm', sprintf('_ch%d.jrc', iChan));
            fprintf('Loading chan %d(%d/%d) from %s\n\t', iChan, iChan1, numel(viChan), P.vcFile_prm);
            t1 = tic;
            vn_ = load_bin_chan_(P, iChan);
            fprintf('\n\ttook %0.1fs\n', toc(t1));
            write_bin_(vcFile_out, vn_);
        end %for
    end
end %func
