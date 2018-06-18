%--------------------------------------------------------------------------
function vl = matchFileExt_(csFiles, vcExt, vlDir)
    % vcExt can be a cell
    % ignore dir
    % matchFileExt_(csFiles, vcExt, vlDir)
    % matchFileExt_(csFiles, csExt, vlDir) %multiple extension check
    if ischar(csFiles), csFiles={csFiles}; end
    vl = false(size(csFiles));

    for i=1:numel(csFiles)
        [~,~,vcExt1] = fileparts(csFiles{i});
        vl(i) = any(strcmpi(vcExt1, vcExt));
    end
    if nargin >= 3
        vl = vl | vlDir; %matches if it's directory
    end
end %func
