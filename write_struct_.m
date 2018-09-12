%--------------------------------------------------------------------------
function write_struct_(vcFile, S)
    % Write a struct S to file vcFile
    try
        warning off
        t1=tic;
        %     S = struct_remove_handles(S); %remove figure handle
        save(vcFile, '-struct', 'S');
        fprintf('Wrote to %s, took %0.1fs.\n', vcFile, toc(t1));
    catch
        fprintf(2, 'Writing struct to file %s failed.\n', vcFile);
    end
end %func
