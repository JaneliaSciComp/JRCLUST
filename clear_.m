%--------------------------------------------------------------------------
% 12/21/17 JJJ: clearing a batch file (v3.2.0)
% 10/15/17 JJJ: clear function memory
% 7/31/17 JJJ: Documentation and test
function clear_(vcFile_prm)
    % Clear JRCLUST global variables
    if nargin<1, vcFile_prm = ''; end
    global fDebug_ui;

    % clear jrc3
    clear(mfilename()); % clear persistent variables in the current file. Same as clear jrc3
    clear global tnWav_spk tnWav_raw trFet_spk mnWav1 mrWav1 mnWav S_gt vcFile_prm_
    clear functions % clear function memory 10/15/17 JJJ

    set(0, 'UserData', []);
    try gpuDevice(1); catch, fprintf(2, 'GPU reset error.\n'); end
    disp('Memory cleared on CPU and GPU');
    fDebug_ui = [];

    if isempty(vcFile_prm), return; end

    % clear specific parameter file. erase prm file related files
    if matchFileExt_(vcFile_prm, '.batch')
        csFile_prm = load_batch_(vcFile_prm);
    elseif matchFileExt_(vcFile_prm, '.prm')
        csFile_prm = {vcFile_prm};
    else
        return;
    end
    for iFile = 1:numel(csFile_prm)
        csFiles_del = strrep(csFile_prm{iFile}, '.prm', {'_jrc.mat', '_spkraw.jrc', '_spkwav.jrc', '_spkfet.jrc', '_log.mat', '_gt1.mat'});
        delete_files_(csFiles_del);
    end
end %func
