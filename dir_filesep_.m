%--------------------------------------------------------------------------
% 11/5/17 JJJ: Created
function vc = dir_filesep_(vc)
    % replace the file seperaation characters
    if isempty(vc), return; end
    vl = vc == '\' | vc == '/';
    if any(vl), vc(vl) = filesep(); end
end %func
