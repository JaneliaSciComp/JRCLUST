function purge(varargin)
    %PURGE Clear global variables and other persistent storage
    % This will eventually be deprecated, but needed while we're transitioning
    global fDebug_ui;

     % clear persistent variables in jrc.m
    clear(fullfile(jrclust.utils.basedir, 'jrc.m'));
    clear global tnWav_spk tnWav_raw trFet_spk mnWav1 mrWav1 mnWav S_gt vcFile_prm_
    clear functions % clear function memory

    set(0, 'UserData', []);
    try
        gpuDevice(1);
    catch
    end

    disp('Memory cleared on CPU and GPU');
    fDebug_ui = [];

    if nargin < 1
        return;
    end

    % clear parameter file (erase prm file related files)
    if matchFileExt_(varargin{1}, '.batch')
        csFile_prm = load_batch_(varargin{1});
    elseif matchFileExt_(varargin{1}, '.prm')
        csFile_prm = {varargin{1}};
    else
        return;
    end
    for iFile = 1:numel(csFile_prm)
        csFiles_del = strrep(csFile_prm{iFile}, '.prm', {'_jrc.mat', '_spkraw.jrc', '_spkwav.jrc', '_spkfet.jrc', '_log.mat', '_gt1.mat'});
        delete_files_(csFiles_del);
    end
end

