%--------------------------------------------------------------------------
% 12/21/17 JJJ: Get the tag by name which is cached (like hash table)
function hObj = get_tag_(vcTag, vcType)
    % clear before starting manual
    % Return from persistent cache
    % Create a new figure if Tag doesn't exist

    persistent S_tag_cache_
    if nargin<2, vcType = []; end
    if isempty(S_tag_cache_)
        S_tag_cache_ = struct();
    else
        if isfield(S_tag_cache_, vcTag)
            hObj = S_tag_cache_.(vcTag);
            if isvalid_(hObj), return; end
        end
    end
    hObj = findobj('Tag', vcTag, 'Type', vcType);
    S_tag_cache_.(vcTag) = hObj;
end %func
