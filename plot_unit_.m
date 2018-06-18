%--------------------------------------------------------------------------
function plot_unit_(S_clu1, hAx, vcColor0)
    if isempty(S_clu1), return; end
    if nargin<2, hAx = axes_new_('FigWav'); end
    if nargin<3, vcColor0 = [0 0 0]; end
    [S0, P, S_clu] = get0_();
    [~, S_figWav] = get_fig_cache_('FigWav');
    maxAmp = S_figWav.maxAmp;
    % plot individual unit
    nSamples = size(S_clu1.mrWav_clu,1);
    vrX = (1:nSamples)'/nSamples;
    vrX([1,end])=nan; % line break

    if ~isequal(vcColor0, [0 0 0])
        trWav1 = zeros(1,1,0);
    else
        trWav1 = S_clu1.trWav;
    end

    % show example traces
    for iWav = size(trWav1,3):-1:0
        if iWav==0
            mrY1 = S_clu1.mrWav_clu / maxAmp;
            lineWidth=1.5;
            vcColor = vcColor0;
        else
            mrY1 = trWav1(:,:,iWav) / maxAmp;
            lineWidth=.5;
            vcColor = .5*[1,1,1];
        end
        vrX1_site = P.mrSiteXY(S_clu1.viSite, 1) / P.um_per_pix;
        vrY1_site = P.mrSiteXY(S_clu1.viSite, 2) / P.um_per_pix;
        mrY1 = bsxfun(@plus, mrY1, vrY1_site');
        mrX1 = bsxfun(@plus, repmat(vrX, [1, size(mrY1, 2)]), vrX1_site');
        line(mrX1(:), mrY1(:), 'Color', vcColor, 'Parent', hAx, 'LineWidth', lineWidth);
    end
    xlabel(hAx, 'X pos [pix]');
    ylabel(hAx, 'Z pos [pix]');
    grid(hAx, 'on');
    xlim_(hAx, [min(mrX1(:)), max(mrX1(:))]);
    ylim_(hAx, [floor(min(mrY1(:))-1), ceil(max(mrY1(:))+1)]);
end %func
