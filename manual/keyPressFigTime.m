%--------------------------------------------------------------------------
function keyPressFigTime(hObject, event, S0)

    if nargin < 3
        S0 = get(0, 'UserData');
    end
    [P, S_clu, hFig] = deal(S0.P, S0.S_clu, hObject);
    figData = get(hFig, 'UserData');

    nSites = numel(P.chanMap);

    switch lower(event.Key)
        case {'leftarrow', 'rightarrow'}
            if ~isVisible_(figData.hAx)
                msgbox_('Channel switching is disabled in the position view');
                return;
            end

            factor = keyModifier(event, 'shift')*3 + 1;
            if strcmpi(event.Key, 'rightarrow')
                figData.primarySite = min(figData.primarySite + factor, nSites);
            else
                figData.primarySite = max(figData.primarySite - factor, 1);
            end

            set(hFig, 'UserData', figData);
            updateFigTime();

        case {'uparrow', 'downarrow'} %change ampl
            if ~isVisible_(figData.hAx)
                msgbox_('Zoom is disabled in the position view'); return;
            end
            rescale_FigTime_(event, S0, P);

        case 'r' %reset view
            if ~isVisible_(figData.hAx)
                return;
            end
            axis_(figData.hAx, [figData.timeLimitSecs, figData.vpp_lim]);
            imrect_set_(figData.hRect, figData.timeLimitSecs, figData.vpp_lim);

        case 'm' %merge
            manualMerge(S0);

        case 'h' %help
            msgbox_(figData.csHelp, 1);

        case 'b' %background spike toggle
            if isVisible_(figData.hAx)
                figData.fPlot0 = toggleVisible_(figData.hPlotBG);
            else
                figData.fPlot0 = toggleVisible_(figData.hPlot0_track);
            end
            set(hFig, 'UserData', figData);

        case 't'
            plotFigTime(S0);

            %     case 'z' % track depth
            %         disp('FigTime:''z'' not implemented yet');
            %         plot_SpikePos_(S0, event);

        case 's' %split. draw a polygon
            if ~isempty(S0.secondarySelectedCluster)
                msgbox_('Select one cluster'); return;
            end
            try
                hPoly = impoly_();
                if isempty(hPoly); return ;end
                mrPolyPos = getPosition(hPoly);
                vrX1 = double(get(figData.hPlotFG, 'XData'));
                vrY1 = double(get(figData.hPlotFG, 'YData'));
                vlIn = inpolygon(vrX1, vrY1, mrPolyPos(:,1), mrPolyPos(:,2));
                hSplit = line(vrX1(vlIn), vrY1(vlIn), 'Color', [1 0 0], 'Marker', '.', 'LineStyle', 'none');
                if strcmpi(userDialog('Split?', 'Confirmation', 'Yes'), 'yes')
                    splitCluster(S0.primarySelectedCluster, vlIn);
                end
                deleteMany(hPoly, hSplit);
            catch
                disp(lasterror());
            end

        case 'p' % update PC projection
            if strcmpi(P.displayFeature, 'kilosort')
                S0.pcTime = mod(S0.pcTime, 3) + 1; % 1 => 2, 2 => 3, 3 => 1
                set(0, 'UserData', S0);
                updateFigTime();
            end


        case 'f' % feature display instead of amplitude display
            if getOr(P, 'fImportKilosort')
                if strcmpi(P.displayFeature, 'vpp')
                    P.displayFeature = 'kilosort';
                else
                    P.displayFeature = 'vpp';
                    set(hFig, 'UserData', figData);
                end
                S0.P = P;
                set(0, 'UserData', S0);

                updateFigTime();
                % also update FigProj
                plotFigProj(S0);
            else
                disp('keyPressFigProj: ''f'': not implemented yet');
            end

        case 'c' % compare pca across channels
            disp('FigTime: Not implemented yet'); return;
            %         hMsg = msgbox_('Plotting...');
            %         figure; hold on;
            %         [mrWav_mean1, viSite1] = mrWav_int_mean_clu_(S0.primarySelectedCluster);
            %         [~, mrPv1] = pca(mrWav_mean1, 'NumComponents', P.nPc_dip, 'Center', 1);
            %         mrPv1 = norm_mr_(mrPv1);
            %
            %         if keyModifier(event, 'control') %show chain of clusters
            %             trPv1 = mrPv1;
            %             iClu_next = get_next_clu_(S_clu, S0.primarySelectedCluster);
            %             viClu_track = S0.primarySelectedCluster;
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
            %         elseif ~isempty(S0.secondarySelectedCluster)
            %             [mrWav_mean2, viSite1] = mrWav_int_mean_clu_(S0.secondarySelectedCluster);
            %             [~, mrPv2] = pca(mrWav_mean2, 'NumComponents', P.nPc_dip);
            %             mrPv2 = match_mrPv_(mrPv2, mrPv1);
            % %             mrPv2 = flip_prinvec_(mrPv2, mrPv1);
            %             mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'k');
            %             mr2plot(norm_mr_(mrPv2), 'scale', 1, 'LineStyle', 'r--');
            %             vcTitle = sprintf('PCA across chan: Clu %d vs %d', S0.primarySelectedCluster, S0.secondarySelectedCluster);
            %         else
            %             mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'r');
            %             vcTitle = sprintf('PCA across chan: Clu %d', S0.primarySelectedCluster);
            %         end
            % %         mr2plot(mrPv1, 'scale', 1, 'LineStyle', 'k');
            %         grid on;
            %         title_(vcTitle);
            % %         if ~isempty(S0.secondarySelectedCluster)
            % %             compare_interp_(Sclu, S0.primarySelectedCluster, S0.secondarySelectedCluster);
            % %         end
            %         try close(hMsg); catch; end

            %     case 'f' %feature export
            %         eval(sprintf('mrFet_clu%d = getFet_clu_(S0.primarySelectedCluster);', S0.primarySelectedCluster));
            %         mrDist1 = squareform(pdist(mrFet1'));
            %         vrFet1 = sqrt(sum(mrFet1.^2));
            %         mrDist1 = bsxfun(@rdivide, mrDist1, vrFet1); %norm
            %         eval(sprintf('assignWorkspace_(mrFet_clu%d);', S0.primarySelectedCluster));

        case 'e' %export selected to workspace
                disp('FigTime: ''e'' not implemented yet'); return;
    end %switch
    % drawnow;
    % figure_wait(0);
end % function
