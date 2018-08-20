%--------------------------------------------------------------------------
function [vrX, vrY, viPlot, tr_dim] = amp2proj_(mrMin, mrMax, maxAmp, maxPair, P)
    if nargin<4, maxPair = []; end
    if nargin<5, P = get0_('P'); end
    % switch lower(P.displayFeature)
    %     case {'vpp', 'vmin', 'vmax'}
    %         mrMax = linmap_(mrMax', [0, maxAmp/2], [0,1], 1);
    %         mrMin = linmap_(mrMin', [0, maxAmp], [0,1], 1);
    %     otherwise
    mrMax = linmap_(mrMax', [0, 1] * maxAmp, [0,1], 1);
    mrMin = linmap_(mrMin', [0, 1] * maxAmp, [0,1], 1);
    % end
    [nEvt, nChans] = size(mrMin);
    if isempty(maxPair), maxPair = nChans; end
    [trX, trY] = deal(nan([nEvt, nChans, nChans], 'single'));
    for chY = 1:nChans
        vrY1 = mrMin(:,chY);
        vlY1 = vrY1>0 & vrY1<1;
        for chX = 1:nChans
            if abs(chX-chY) > maxPair, continue; end
            if chY > chX
                vrX1 = mrMin(:,chX);
            else
                vrX1 = mrMax(:,chX);
            end
            viPlot1 = find(vrX1>0 & vrX1<1 & vlY1);
            trX(viPlot1,chY,chX) = vrX1(viPlot1) + chX - 1;
            trY(viPlot1,chY,chX) = vrY1(viPlot1) + chY - 1;
        end
    end
    % plot projection
    viPlot = find(~isnan(trX) & ~isnan(trY));
    vrX = trX(viPlot);  vrX=vrX(:);
    vrY = trY(viPlot);  vrY=vrY(:);
    tr_dim = size(trX);
end %func
