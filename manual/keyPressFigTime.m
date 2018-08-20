%--------------------------------------------------------------------------
function keyPressFigTime(hObject, event, S0)

    if nargin<3, S0 = get(0, 'UserData'); end
    [P, S_clu, hFig] = deal(S0.P, S0.S_clu, hObject);
    S_fig = get(hFig, 'UserData');

    nSites = numel(P.chanMap);
    % set(hObject, 'Pointer', 'watch');
    % figure_wait(1);
    switch lower(event.Key)
        case {'leftarrow', 'rightarrow'}
            if ~isVisible_(S_fig.hAx)
                msgbox_('Channel switching is disabled in the position view'); return;
            end
            factor = keyModifier(event, 'shift')*3 + 1;
            if strcmpi(event.Key, 'rightarrow')
                S_fig.iSite = min(S_fig.iSite + factor, nSites);
            else
                S_fig.iSite = max(S_fig.iSite - factor, 1);
            end
            set(hFig, 'UserData', S_fig);
            updateFigTime();

        case {'uparrow', 'downarrow'} %change ampl
            if ~isVisible_(S_fig.hAx)
                msgbox_('Zoom is disabled in the position view'); return;
            end
            rescale_FigTime_(event, S0, P);

        case 'r' %reset view
            if ~isVisible_(S_fig.hAx), return; end
            axis_(S_fig.hAx, [S_fig.time_lim, S_fig.vpp_lim]);
            imrect_set_(S_fig.hRect, S_fig.time_lim, S_fig.vpp_lim);

        case 'm' %merge
            ui_merge_(S0);

        case 'h' %help
            msgbox_(S_fig.csHelp, 1);

        case 'b' %background spike toggle
                if isVisible_(S_fig.hAx)
                    S_fig.fPlot0 = toggleVisible_(S_fig.hPlot0);
                else
                    S_fig.fPlot0 = toggleVisible_(S_fig.hPlot0_track);
                end
                set(hFig, 'UserData', S_fig);

        case 't'
                plotFigTime(S0);

                %     case 'z' % track depth
                %         disp('FigTime:''z'' not implemented yet');
                %         plot_SpikePos_(S0, event);

        case 's' %split. draw a polygon
                if ~isempty(S0.iCluPaste)
                    msgbox_('Select one cluster'); return;
                end
                try
                    hPoly = impoly_();
                    if isempty(hPoly); return ;end
                    mrPolyPos = getPosition(hPoly);
                    vrX1 = double(get(S_fig.hPlot1, 'XData'));
                    vrY1 = double(get(S_fig.hPlot1, 'YData'));
                    vlIn = inpolygon(vrX1, vrY1, mrPolyPos(:,1), mrPolyPos(:,2));
                    hSplit = line(vrX1(vlIn), vrY1(vlIn), 'Color', [1 0 0], 'Marker', '.', 'LineStyle', 'none');
                    if strcmpi(userDialog('Split?', 'Confirmation', 'Yes'), 'yes')
                        split_clu_(S0.iCluCopy, vlIn);
                    end
                    deleteMany(hPoly, hSplit);
                catch
                    disp(lasterror());
                end

        case 'p' % update projection view
                vrPos = getPosition(S_fig.hRect);
                tlim_proj = [vrPos(1), sum(vrPos([1,3]))];
                P.tlim_proj = tlim_proj;
                plotFigProj(S0);

        case 'k' % kilosort amplitude scaling factor/spike
            if getOr(P, 'fImportKilosort', 0)
                displayFeature = getOr(P, 'displayFeature', 'vpp');
                if strcmpi(displayFeature, 'vpp')
                    P.displayFeature = 'kilosort';
                elseif strcmpi(displayFeature, 'kilosort')
                    P.displayFeature = 'vpp'; % others not implemented yet
                end

                setUserData(P);
                updateFigTime();
            end

                %     case 'f' % feature display instead of amplitude display
                %         if strcmpi(P.displayFeature, 'fet')
                %             P.displayFeature = 'vpp';
                %         else
                %             P.displayFeature = 'fet';
                %         end
                %         setUserData(P);
                %         updateFigTime();

        case 'c' % compare pca across channels
                disp('FigTime: Not implemented yet'); return;
                %         hMsg = msgbox_('Plotting...');
                %         figure; hold on;
                %         [mrWav_mean1, viSite1] = mrWav_int_mean_clu_(S0.iCluCopy);
                %         [~, mrPv1] = pca(mrWav_mean1, 'NumComponents', P.nPc_dip, 'Center', 1);
                %         mrPv1 = norm_mr_(mrPv1);
                %
                %         if keyModifier(event, 'control') %show chain of clusters
                %             trPv1 = mrPv1;
                %             iClu_next = get_next_clu_(S_clu, S0.iCluCopy);
                %             viClu_track = S0.iCluCopy;
                %             while ~isempty(iClu_next)
                %                 [mrWav_mean1, viSite1] = mrWav_int_mean_clu_(iClu_next);
                %                 [~, mrPv1a] = pca(mrWav_mean1, 'NumComponents', P.nPc_dip, 'Center', 1);
                %                 mrPv1a = norm_mr_(mrPv1a);
                %                 mrPv1 = flip_prinvec_(mrPv1a, mean(trPv1,3));
                %                 trPv1 = cat(3, trPv1, mrPv1);
                %                 viClu_track(end+1) = iClu_next;
                %
                %                 iClu_next = get_next_clu_(S_clu, iClu_next);
                %             end
                %             multiplot(plot(nan,nan,'k'), 1, 1:size(trPv1,1), trPv1);
                % %             mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'k');
                %             vcTitle = sprintf('PCA across chan: Clu %s', sprintf('%d,', viClu_track));
                %         elseif ~isempty(S0.iCluPaste)
                %             [mrWav_mean2, viSite1] = mrWav_int_mean_clu_(S0.iCluPaste);
                %             [~, mrPv2] = pca(mrWav_mean2, 'NumComponents', P.nPc_dip);
                %             mrPv2 = match_mrPv_(mrPv2, mrPv1);
                % %             mrPv2 = flip_prinvec_(mrPv2, mrPv1);
                %             mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'k');
                %             mr2plot(norm_mr_(mrPv2), 'scale', 1, 'LineStyle', 'r--');
                %             vcTitle = sprintf('PCA across chan: Clu %d vs %d', S0.iCluCopy, S0.iCluPaste);
                %         else
                %             mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'r');
                %             vcTitle = sprintf('PCA across chan: Clu %d', S0.iCluCopy);
                %         end
                % %         mr2plot(mrPv1, 'scale', 1, 'LineStyle', 'k');
                %         grid on;
                %         title_(vcTitle);
                % %         if ~isempty(S0.iCluPaste)
                % %             compare_interp_(Sclu, S0.iCluCopy, S0.iCluPaste);
                % %         end
                %         try close(hMsg); catch; end

                %     case 'f' %feature export
                %         eval(sprintf('mrFet_clu%d = getFet_clu_(S0.iCluCopy);', S0.iCluCopy));
                %         mrDist1 = squareform(pdist(mrFet1'));
                %         vrFet1 = sqrt(sum(mrFet1.^2));
                %         mrDist1 = bsxfun(@rdivide, mrDist1, vrFet1); %norm
                %         eval(sprintf('assignWorkspace_(mrFet_clu%d);', S0.iCluCopy));

        case 'e' %export selected to workspace
                disp('FigTime: ''e'' not implemented yet'); return;
    end %switch
    % drawnow;
    % figure_wait(0);
end %func
