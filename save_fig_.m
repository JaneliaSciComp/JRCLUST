%--------------------------------------------------------------------------
function save_fig_(vcFile_png, hFig, fClose);
    % Save a figure to a figure file
    if nargin<2, hFig = []; end
    if nargin<3, fClose = 0; end
    if isempty(hFig), hFig = gcf; end

    % vcFile_png = strrep(P.prmFile, '.prm', '.png');
    try
        hMsg = msgbox_('Saving figure... (this closes automaticall)');
        drawnow;
        saveas(hFig, vcFile_png);
        fprintf('Saved figure to %s.\n', vcFile_png);
        if fClose, tryClose(hFig); end
        tryClose(hMsg);
    catch
        fprintf(2, 'Failed to save a figure to %s.\n', vcFile_png);
    end
end %func
