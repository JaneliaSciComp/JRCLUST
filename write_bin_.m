%--------------------------------------------------------------------------
function dimm_mr = write_bin_(vcFile, mr)
    t1=tic;
    dimm_mr = size(mr);
    fVerbose = 1;
    if isempty(mr), return; end
    if isstruct(mr)
        save(vcFile, '-struct', 'mr', '-v7.3'); % save to matlab file
    else
        if ischar(vcFile)
            fid_w = fopen(vcFile, 'W');
        else
            fid_w = vcFile;
        end
        fwrite(fid_w, mr, class(mr));
        if ischar(vcFile)
            fclose(fid_w);
        else
            fVerbose = 0;
        end
    end
    if fVerbose
        fprintf('Writing to %s took %0.1fs\n', vcFile, toc(t1));
    end
end %func
