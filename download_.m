%--------------------------------------------------------------------------
function download_(vcMode)
    S_cfg = read_cfg_();
    % if strcmpi(pwd(), S_cfg.path_alpha), disp('cannot overwrite alpha'); return; end
    switch lower(vcMode)
        case {'sample', 'neuropix2', 'neuropixels2', 'phase2', 'phaseii'}
        csLink = S_cfg.path_sample_phase2;
        case {'sample3', 'neuropix3' 'neuropixels3', 'phase3', 'phaseiii'}
        csLink = S_cfg.path_sample_phase3;
        otherwise
        disp('Invalid selection. Try "jrc download sample or jrc download sample3".');
        return;
    end %switch

    t1 = tic;
    fprintf('Downloading sample files. This can take up to several minutes.\n');
    vlSuccess = download_files_(csLink);
    fprintf('\t%d/%d files downloaded. Took %0.1fs\n', ...
    sum(vlSuccess), numel(vlSuccess), toc(t1));
end %func
