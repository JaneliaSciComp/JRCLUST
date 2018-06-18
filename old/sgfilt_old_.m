%--------------------------------------------------------------------------
% 12/21/17 JJJ: Annotate on the x-axis intead of floating text box (faster)
% function vhText = text_nClu_(S_clu, hAx)
% % update x label
%
% nClu = numel(S_clu.vnSpk_clu);
% viSite_clu = double(S_clu.viSite_clu);
% y_offset = .3;
% vhText = zeros(nClu+1,1);
% for iClu=1:nClu
%     n1 = S_clu.vnSpk_clu(iClu);
% %         vcText1 = sprintf('%d, %0.1f', n1, n1/t_dur);
%     vcText1 = sprintf('%d', n1); %show numbers
%     vhText(iClu) = text(iClu, viSite_clu(iClu) + y_offset, vcText1, ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Parent', hAx);
% end
% vhText(end) = text(0,0,'#spk', 'VerticalAlignment', 'bottom', 'Parent', hAx);
% end %func


%--------------------------------------------------------------------------
function mnWav2 = sgfilt_old_(mnWav, nDiff_filt, fGpu)
    % works for a vector, matrix and tensor

    if nargin<3, fGpu = isGpu_(mnWav); end
    n1 = size(mnWav,1);
    if n1==1, n1 = size(mnWav,2);  end
    [via1, via2, via3, via4, vib1, vib2, vib3, vib4] = sgfilt4_(n1, fGpu);

    if isvector(mnWav)
        switch nDiff_filt
            case 1
            mnWav2 = mnWav(via1) - mnWav(vib1);
            case 2
            mnWav2 = 2*(mnWav(via2) - mnWav(vib2)) + mnWav(via1) - mnWav(vib1);
            case 3
            mnWav2 = 3*(mnWav(via3) - mnWav(vib3)) + 2*(mnWav(via2) - mnWav(vib2)) + mnWav(via1) - mnWav(vib1);
            otherwise
            mnWav2 = 4*(mnWav(via4) - mnWav(vib4)) + 3*(mnWav(via3) - mnWav(vib3)) + 2*(mnWav(via2) - mnWav(vib2)) + mnWav(via1) - mnWav(vib1);
        end %switch
    elseif ismatrix(mnWav)
        switch nDiff_filt
            case 1
            mnWav2 = mnWav(via1,:) - mnWav(vib1,:);
            case 2
            mnWav2 = 2*(mnWav(via2,:) - mnWav(vib2,:)) + mnWav(via1,:) - mnWav(vib1,:);
            case 3
            mnWav2 = 3*(mnWav(via3,:) - mnWav(vib3,:)) + 2*(mnWav(via2,:) - mnWav(vib2,:)) + mnWav(via1,:) - mnWav(vib1,:);
            otherwise
            mnWav2 = 4*(mnWav(via4,:) - mnWav(vib4,:)) + 3*(mnWav(via3,:) - mnWav(vib3,:)) + 2*(mnWav(via2,:) - mnWav(vib2,:)) + mnWav(via1,:) - mnWav(vib1,:);
        end %switch
    else
        switch nDiff_filt
            case 1
            mnWav2 = mnWav(via1,:,:) - mnWav(vib1,:,:);
            case 2
            mnWav2 = 2*(mnWav(via2,:,:) - mnWav(vib2,:,:)) + mnWav(via1,:,:) - mnWav(vib1,:,:);
            case 3
            mnWav2 = 3*(mnWav(via3,:,:) - mnWav(vib3,:,:)) + 2*(mnWav(via2,:,:) - mnWav(vib2,:,:)) + mnWav(via1,:,:) - mnWav(vib1,:,:);
            otherwise
            mnWav2 = 4*(mnWav(via4,:,:) - mnWav(vib4,:,:)) + 3*(mnWav(via3,:,:) - mnWav(vib3,:,:)) + 2*(mnWav(via2,:,:) - mnWav(vib2,:,:)) + mnWav(via1,:,:) - mnWav(vib1,:,:);
        end %switch
    end
end %func
