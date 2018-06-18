%--------------------------------------------------------------------------
% 12/5/17 JJJ: Colors are hard coded
% 9/17/17 JJJ: nGroups fixed, can have less than 7 clusters
function vhPlot = plot_group_(hAx, mrX, mrY, varargin)
    mrColor = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840]';
    nGroups = min(size(mrColor,2), size(mrX,2));
    mrColor = mrColor(:,1:nGroups);
    vhPlot = zeros(nGroups, 1);
    hold(hAx, 'on');
    for iGroup=1:nGroups
        vrX1 = mrX(:, iGroup:nGroups:end);
        vrY1 = mrY(:, iGroup:nGroups:end);
        vhPlot(iGroup) = plot(hAx, vrX1(:), vrY1(:), varargin{:}, 'Color', mrColor(:,iGroup)');
    end
end %func
