%--------------------------------------------------------------------------
function fExit = save_manual_(varargin)
    % TODO: remove figure handles from S0
    if nargin==1
        P = varargin{1};
    else
        P = get0_('P');
    end
    vcFile_jrc = subsFileExt_(P.vcFile_prm, '_jrc.mat');
    fExit = 1;
    switch lower(questdlg_(['Save to ', vcFile_jrc, ' ?'], 'Confirmation', 'Yes'))
        case 'yes'
        hMsg = jrclust.utils.qMsgBox('Saving... (this closes automatically)');
        save0_(vcFile_jrc); % 1 will skip figure saving
        fExit = 1;
        close_(hMsg);
        case 'no'
        fExit = 1;
        case 'cancel'
        fExit = 0;
        return;
    end
end %func;
