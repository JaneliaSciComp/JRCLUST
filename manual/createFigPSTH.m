%--------------------------------------------------------------------------
function [hFig, hFig_b] = createFigPSTH(hFig, hFig_b, P, nStims)

    % Figure handle for the primarySelectedCluster
    [axoffset, axlen] = deal(.08, 1/nStims);

    if ~tryIsValid(hFig)
        hFig = createFigure('FigTrial', [.5  .5 .5 .5], P.trialFile, 0, 0);
        [vhAx1, vhAx2] = deal(nan(nStims, 1));
        for iStim = 1:nStims
            axoffset_ = axoffset + (iStim-1) * axlen;
            vhAx1(iStim) = axes('Parent', hFig, 'Position',[.08 axoffset_ .9 axlen*.68]);
            vhAx2(iStim) = axes('Parent', hFig, 'Position',[.08 axoffset_ + axlen*.68 .9 axlen*.2]);
        end
        vcColor = 'k';
        set(hFig, 'UserData', makeStruct_(vhAx1, vhAx2, vcColor));
    end

    % Figure handle for the secondarySelectedCluster
    if ~tryIsValid(hFig_b)
        hFig_b = createFigure('FigTrial_b', [.5  0 .5 .5], P.trialFile, 0, 0);
        set(hFig_b, 'Visible', 'off');
        [vhAx1, vhAx2] = deal(nan(nStims, 1));
        for iStim = 1:nStims
            axoffset_ = axoffset + (iStim-1) * axlen;
            vhAx1(iStim) = axes('Parent', hFig_b, 'Position',[.08 axoffset_ .9 axlen*.68]);
            vhAx2(iStim) = axes('Parent', hFig_b, 'Position',[.08 axoffset_ + axlen*.68 .9 axlen*.2]);
        end
        vcColor = 'r';
        set(hFig_b, 'UserData', makeStruct_(vhAx1, vhAx2, vcColor));
    end
end % function
