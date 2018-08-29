%--------------------------------------------------------------------------
% function P = appendStruct_(P, varargin)
% % backward compatibility
% P = mergeStructs(P, varargin{:});
% end


%--------------------------------------------------------------------------
function vcPath = replacePath_(vcPath1, vcPath2)
    % replace path1 with path2
    [~, vcFname1, vcExt1] = fileparts(vcPath1);
    [vcDir2,~,~] = fileparts(vcPath2);
    if ~isempty(vcDir2)
        vcPath = [vcDir2, filesep(), vcFname1, vcExt1];
    else
        vcPath = [vcFname1, vcExt1];
    end
end % function
