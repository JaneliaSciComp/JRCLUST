%--------------------------------------------------------------------------
function setpath_()
    % Reset to the Matlab default path in case user overrided the system default functions
    persistent fPathSet
    if fPathSet, return; end

    % restoredefaultpath;
    % disp('Matlab path is temporarily reset to the factory default for this session.');

    % [vcPath_jrc, ~, ~] = fileparts(mfilename('fullpath')); % get current directory
    % addpath(vcPath_jrc);
    fPathSet = 1;
end
