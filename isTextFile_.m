%--------------------------------------------------------------------------
% 8/4/17 JJJ: Function created, tested and documented
function flag = isTextFile_(vcFile_txt, textFileExtensions)
    % Return true if a file is a text file
    % Accept list of text files (textFileExtensions}, provided by extensions

    if nargin < 2
        textFileExtensions = {'.txt', '.batch'};
    end

    if ischar(vcFile_txt)
        flag = matchFileExt_(vcFile_txt, textFileExtensions);
        %     flag = ~isempty(regexp(lower(vcFile_txt), '.txt$'));
    else
        flag = 0;
    end
end %func
