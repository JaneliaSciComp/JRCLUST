%--------------------------------------------------------------------------
function export_quality_(varargin)
    % export_csv_(hObject, event)
    if nargin==1
        P = varargin{1};
        fGui = 0;
    else
        P = get0_('P');
        fGui = 1;
    end

    [S0, P] = load_cached_(P);
    if ~isfield(S0, 'S_clu'), fprintf(2, 'File must be sorted first.\n'); return; end
    S = S0.S_clu;
    [unit_id, SNR, center_site, nSpikes, xpos, ypos, uV_min, uV_pp, IsoDist, LRatio, IsiRat, note] = ...
    deal((1:S.nClu)', S.vrSnr_clu(:), S.clusterSites(:), S.nSpikesPerCluster(:), ...
    S.clusterXPositions(:), S.clusterXPositions(:), S.vrVmin_uv_clu, S.vrVpp_uv_clu, ...
    S.vrIsoDist_clu(:), S.vrLRatio_clu(:), S.vrIsiRatio_clu(:), S.clusterNotes(:));
    %note(cellfun(@isempty, note)) = '';
    table_ = table(unit_id, SNR, center_site, nSpikes, xpos, ypos, uV_min, uV_pp, IsoDist, LRatio, IsiRat, note);
    disp(table_);

    vcFile_csv = subsFileExt_(P.paramFile, '_quality.csv');
    writetable(table_, vcFile_csv);
    csMsg = { ...
    sprintf('Wrote to %s. Columns:', vcFile_csv), ...
    sprintf('\tColumn 1: unit_id: Unit ID'), ...
    sprintf('\tColumn 2: SNR: |Vp/Vrms|; Vp: negative peak amplitude of the peak site; Vrms: SD of the Gaussian noise (estimated from MAD)'), ...
    sprintf('\tColumn 3: center_site: Peak site number which contains the most negative peak amplitude'), ...
    sprintf('\tColumn 4: nSpikes: Number of spikes'), ...
    sprintf('\tColumn 5: xpos: x position (width dimension) center-of-mass'), ...
    sprintf('\tColumn 6: ypos: y position (depth dimension) center-of-mass, referenced from the tip'), ...
    sprintf('\tColumn 7: uV_min: negative peak voltage (microvolts)'), ...
    sprintf('\tColumn 8: uV_pp: peak-to-peak voltage (microvolts)'), ...
    sprintf('\tColumn 9: IsoDist: Isolation distance quality metric'), ...
    sprintf('\tColumn 10: LRatio: L-ratio quality metric'), ...
    sprintf('\tColumn 11: IsiRat: ISI-ratio quality metric'), ...
    sprintf('\tColumn 12: note: user comments')};

    cellfun(@(x)fprintf('%s\n',x), csMsg);
    if fGui, msgbox_(csMsg); end
end %func
