%--------------------------------------------------------------------------
function P = file2struct__(vcFile_file2struct)
    % James Jun 2017 May 23
    % Run a text file as .m script and result saved to a struct P
    % _prm and _prb can now be called .prm and .prb files

    % load text file. trim and line break. remove comments.  replace
    csLines_file2struct = file2lines_(vcFile_file2struct);
    csLines_file2struct = strip_comments_(csLines_file2struct);
    if isempty(csLines_file2struct), P=[]; return; end

    try
        eval(cell2mat(csLines_file2struct'));

        S_ws = whos();
        csVars = {S_ws.name};
        csVars = setdiff(csVars, {'csLines_file2struct', 'vcFile_file2struct'});
        for i=1:numel(csVars)
            eval(sprintf('a = %s;', csVars{i}));
            P.(csVars{i}) = a;
        end
    catch
        disp(lasterr());
        P=[];
    end
end % function
