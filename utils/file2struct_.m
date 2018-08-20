%--------------------------------------------------------------------------
function P = file2struct_(vcFile_file2struct)
    % Run a text file as .m script and result saved to a struct P
    % _prm and _prb can now be called .prm and .prb files

    try
        P = file2struct__(vcFile_file2struct); % new version
    catch
        P = file2struct_1_(vcFile_file2struct); % old version
    end
    % if isempty(P)
    %     copyfile(vcFile_file2struct, 'temp_eval.m', 'f');
    %     try
    %         eval('temp_eval.m');
    %     catch
    %         disp(lasterr());
    %     end
    % end
end %func
