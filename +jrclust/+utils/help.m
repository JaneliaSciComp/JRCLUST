function help(varargin)
    %HELP Display the documentation
    md = jrclust.utils.info();

    try
        web(md.docsSite);
    catch ME
        fprintf('Please visit %s to view the documentation\n', md.docsSite);
    end
end
