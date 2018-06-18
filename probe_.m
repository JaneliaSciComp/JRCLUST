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
    % S_prb = file2struct_(vcFile_prb);
    S_prb = load_prb_(vcFile_prb);
    % if ~isfield(S_prb, 'shank'), S_prb.shank = ones(size(S_prb.channels)); end

    % hFig = figure; hold on;
    hFig = create_figure_('FigProbe', [0 0 .5 1], vcFile_prb);
    % viShank = unique(S_prb.shank);
    % vcColor_shank = 'kbgrcm'; % up to
    % for iShank=1:numel(viShank)
    %     viSite1 = find(S_prb.shank == iShank);
    hPatch = plot_probe_(S_prb.mrSiteXY, S_prb.vrSiteHW, S_prb.viSite2Chan, S_prb.viShank_site);
    %     iShank1 = mod(iShank-1, numel(vcColor_shank))+1;
    %     set(hPatch, 'EdgeColor', vcColor_shank(iShank1));
    % end
    % plot_probe_(mrSiteXY, vrSiteHW, viSite2Chan, vrVpp, hFig)
    % vrPos0 = get(0, 'ScreenSize');
    % set(hFig, 'OuterPosition', vrPos0, 'Color', 'w');
    axis equal;
    % set(hFig, 'Name', vcFile_prb, 'NumberTitle', 'off');
    edit(vcFile_prb); %show probe file
    figure(hFig);
end %func
