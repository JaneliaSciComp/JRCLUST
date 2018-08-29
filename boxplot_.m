%--------------------------------------------------------------------------
function boxplot_(vrY, vrX, xbin, xlim1)
    % range to plot: xlim1
    vcMode = 'both';
    fOdd = 1;
    if fOdd
        viX = ceil((vrX+xbin/2)/xbin);
    else
        viX = ceil(vrX/xbin);
    end
    ilim = ceil(xlim1/xbin);
    viX(viX<ilim(1))=ilim(1);
    viX(viX>ilim(end))=ilim(end);

    nbins = diff(ilim)+1;
    mrYp = zeros(nbins,3);
    viXp = (ilim(1):ilim(end));
    if fOdd
        vrXp = viXp * xbin - xbin;
    else
        vrXp = viXp * xbin - xbin/2;
    end
    for ibin=1:nbins
        try
            ibin1 = viXp(ibin);
            mrYp(ibin, :) = quantile(vrY(viX==ibin1), [.25, .5, .75]);
        catch
            ;
        end
    end

    switch lower(vcMode)
        case 'stairs'
        vrXp = viXp - xbin/2;
        vrXp(end+1)=vrXp(end)+xbin;
        mrYp(end+1,:) = mrYp(end,:);
        stairs(vrXp, mrYp); grid on;
        case 'line'
        plot(vrXp, mrYp); grid on;
        case 'both'
        hold on;
        vrXp1 = vrXp - xbin/2;
        vrXp1(end+1)=vrXp1(end)+xbin;
        stairs(vrXp1, mrYp([1:end,end],[1, 3]), 'k-');

        plot(vrXp, mrYp(:,2), 'k.-', 'LineWidth', 1);
        grid on;
    end
end % function
