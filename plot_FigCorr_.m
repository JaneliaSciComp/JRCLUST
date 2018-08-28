%--------------------------------------------------------------------------
function plot_FigCorr_(S0)
    % hFigCorr plot
    jitter_ms = .5; % bin size for correlation plot
    nLags_ms = 25; %show 25 msec

    if nargin<1, S0 = get(0, 'UserData'); end
    P = S0.P; S_clu = S0.S_clu;
    P.jitter_ms = jitter_ms;
    P.nLags_ms = nLags_ms;

    [hFig, S_fig] = getCachedFig('FigCorr');
    iClu1 = S0.primarySelectedCluster;
    iClu2 = S0.secondarySelectedCluster;
    if isempty(iClu2), iClu2 = iClu1; end

    jitter = round(P.sampleRateHz / 1000 * P.jitter_ms); %0.5 ms
    nLags = round(P.nLags_ms / P.jitter_ms);

    vi1 = int32(double(S_clu_time_(S_clu, iClu1)) /jitter);

    if iClu1~=iClu2
        vi1 = [vi1, vi1-1, vi1+1]; %allow missing one
    end
    vi2 = int32(double(S_clu_time_(S_clu, iClu2)) /jitter);
    viLag = -nLags:nLags;
    vnCnt = zeros(size(viLag));
    for iLag=1:numel(viLag)
        if iClu1 == iClu2 && viLag(iLag)==0, continue; end
        vnCnt(iLag) = numel(intersect(vi1, vi2+viLag(iLag)));
    end
    vrTime_lag = viLag * P.jitter_ms;

    %--------------
    % draw
    if isempty(S_fig)
        S_fig.hAx = newAxes(hFig);
        S_fig.hBar = bar(vrTime_lag, vnCnt, 1);
        xlabel('Time (ms)');
        ylabel('Counts');
        grid on;
        set(S_fig.hAx, 'YScale', 'log');
    else
        set(S_fig.hBar, 'XData', vrTime_lag, 'YData', vnCnt);
    end
    title_(S_fig.hAx, sprintf('Clu%d vs Clu%d', iClu1, iClu2));
    xlim_(S_fig.hAx, [-nLags, nLags] * P.jitter_ms);
    set(hFig, 'UserData', S_fig);
end %func
