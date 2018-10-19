%--------------------------------------------------------------------------
function probe_(vcFile_prb)
    % if nargin<1, vcFile_prb='imec2.prb'; end

    % set prb file
    if matchFileExt_(vcFile_prb, {'.bin', '.dat'})
        vcFile_prb = subsFileExt_(vcFile_prb, '.prm');
    end
    if matchFileExt_(vcFile_prb, '.prm')
        vcFile_prm = vcFile_prb;
        P = loadParam_(vcFile_prm);
        vcFile_prb = P.probe_file;
        if ~exist(vcFile_prb, 'file')
            vcFile_prb = replacePath_(vcFile_prb, vcFile_prm);
        end
    end
    vcFile_prb = find_prb_(vcFile_prb);
    S_prb = load_prb_(vcFile_prb);

    hFig = create_figure_('FigProbe', [0 0 .5 1], vcFile_prb);
    hPatch = plot_probe_(S_prb.mrSiteXY, S_prb.vrSiteHW, S_prb.viSite2Chan, S_prb.viShank_site);
    axis equal;
    edit(vcFile_prb); %show probe file
    figure(hFig);
end %func

%% local functions
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
