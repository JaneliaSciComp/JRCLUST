%--------------------------------------------------------------------------
% convert a cell string to a character array separated by '\n'.
% 17/7/26 code cleanup and testing
function vc = toString_(cs)
    if isempty(cs)
        vc = '';
    elseif ischar(cs)
        vc = cs(:)';
    elseif iscell(cs)
        vc = sprintf('%s\n', cs{:});
    elseif isnumeric(cs)
        vc = sprintf('%f, ', cs(:));
    end
end %func
