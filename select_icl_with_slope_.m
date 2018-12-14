function [ icl ] = select_icl_with_slope_(S_clu, thresh)
    %Calculate regression line to select cluster centers on each site
    % Hidehiko Inagako, 20160425

    if nargin < 2
        thresh = 0.2;
    end

    global fDebug_ui;

    x_range=-4:0.1:1.5;
    y_range=-3.5:0.1:1.5;

    vrLR = log10(S_clu.rho);
    vrLD = log10(S_clu.delta);

    cviSpk_site = get0_('cviSpk_site');
    if isempty(cviSpk_site)
        viSite_spk = get0_('viSite_spk');
        cviSpk_site = arrayfun(@(iSite) find(viSite_spk == iSite), 1:max(viSite_spk), 'UniformOutput', 0);
    end

    nSites = numel(cviSpk_site);
    cviIcl_site = cell(nSites, 1);

    for iSite = 1:nSites    
        viSpk1 = cviSpk_site{iSite};
        if isempty(viSpk1)
            continue;
        end

        vrLD1 = vrLD(viSpk1);
        vrLR1 = vrLR(viSpk1);

        % acquire density of scatter 
        I2 = hist3([vrLR1(:), vrLD1(:)],'edges',{x_range,y_range});
        I2 = medfilt2(I2');

        up_contour_x = x_range;
        up_contour_y = nan(size(x_range));
        for iX = 1:size(I2, 2) % find largest delta bin for each rho bin
            idxLargest = find(I2(:,iX)>1,1,'last');
            if ~isempty(idxLargest)
                up_contour_y(iX) = y_range(idxLargest);
            end
        end

        % 0 or 1 delta bins found for these rho bins
        viKill = isnan(up_contour_y);
        up_contour_x(viKill) = [];
        up_contour_y(viKill) = [];

        if fDebug_ui % plot upper boundary for inspection
            figure; hold on;
            imagesc(I2, 'XData', x_range, 'YData', y_range); 
            axis xy;
            hold on; plot(vrLR1, vrLD1, '.');
            plot(up_contour_x, up_contour_y, 'r-');
            pause;
        end

        % regression of upper edge line
        X = [ones(length(up_contour_x),1) up_contour_x'];
        b = X \ up_contour_y';
        vrLDF1 = b(1) + b(2)*vrLR1 + thresh; % fit line
        cviIcl_site{iSite} = viSpk1(vrLD1 > vrLDF1);

        % get the outliers

    %     x_hist=-4:0.01:2;
    %     %figure;hist(data2(:,2),x_hist)
    %     counts=hist(data2(:,2),x_hist);
    % 
    %     counts2=counts(x_hist>0);
    %     first_0=find(counts2==0,1);
    %     tmp_thr=first_0*0.01;
    %     thresh=max([thresh,tmp_thr]);

        % figure
        % hold on
        % plot(data2(:,1),data2(:,2),'b.')
        % plot(data2((data2(:,2)>thresh),1),data2((data2(:,2)>thresh),2),'r.')
        % hold off


    %     peak_id=id_of_pints(data2(:,2)>thresh);
    %     icl=[icl,peak_id];

    end

    try
        icl = cell2mat(cviIcl_site);
    catch
        icl = cell2mat(cviIcl_site');
    end

    if fDebug_ui % plot centers on all sites for inspection
        figure;
        plot(vrLR(1:10:end), vrLD(1:10:end), '.', vrLR(icl), vrLD(icl), 'r.');
    end
end

