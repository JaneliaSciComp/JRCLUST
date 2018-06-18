%--------------------------------------------------------------------------
function hPatch = plot_probe_(mrSiteXY, vrSiteHW, viSite2Chan, vrVpp, hFig)
    if nargin<3, viSite2Chan=[]; end
    if nargin<4, vrVpp=[]; end
    if nargin<5, hFig=[]; end

    if isempty(hFig)
        hFig = gcf;
    else
        figure(hFig);
    end

    vrX = [0 0 1 1] * vrSiteHW(2);
    vrY = [0 1 1 0] * vrSiteHW(1);

    mrPatchX = bsxfun(@plus, mrSiteXY(:,1)', vrX(:));
    mrPatchY = bsxfun(@plus, mrSiteXY(:,2)', vrY(:));
    nSites = size(mrSiteXY,1);
    if ~isempty(vrVpp)
        hPatch = patch(mrPatchX, mrPatchY, repmat(vrVpp(:)', [4, 1]), 'EdgeColor', 'k', 'FaceColor', 'flat'); %[0 0 0], 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', [0 0 0], 'FaceAlpha', 0);
        caxis([0, max(vrVpp)]);
        colormap jet;
    else
        hPatch = patch(mrPatchX, mrPatchY, 'w', 'EdgeColor', 'k'); %[0 0 0], 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', [0 0 0], 'FaceAlpha', 0);
    end
    if ~isempty(viSite2Chan)
        csText = arrayfun(@(i)sprintf('%d/%d', i, viSite2Chan(i)), 1:numel(viSite2Chan), 'UniformOutput', 0);
    else
        csText = arrayfun(@(i)sprintf('%d', i), 1:nSites, 'UniformOutput', 0);
    end
    hText = text(mrSiteXY(:,1), mrSiteXY(:,2), csText, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    axis_([min(mrPatchX(:)), max(mrPatchX(:)), min(mrPatchY(:)), max(mrPatchY(:))]);
    title('Site# / Chan# (zoom: wheel; pan: hold wheel & drag)');
    % vrPos = get(gcf, 'Position');
    xlabel('X Position (\mum)');
    ylabel('Y Position (\mum)');
    mouse_figure(hFig);
end
