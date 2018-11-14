%--------------------------------------------------------------------------
function mnWav = ndiff_(mnWav, nDiff_filt)
    % works for a vector, matrix and tensor. Like sgdiff but averaging instead
    % fInvert_filter = 0;

    % n1 = size(mnWav,1);
    % if n1==1, n1 = size(mnWav,2);  end
    % if nDiff_filt==0, mnWav1 = mnWav; return; end
    % if fInvert_filter
    %     [miB, miA] = sgfilt_init_(n1, nDiff_filt, isGpu_(mnWav));
    % else
    %     [miA, miB] = sgfilt_init_(n1, nDiff_filt, isGpu_(mnWav));
    % end
    if nDiff_filt==0, return; end
    % try
    switch nDiff_filt
        case 1
        mnWav(1:end-1,:) = diff(mnWav);
        mnWav(end,:) = 0;
        case 2
        mnWav(2:end-2,:) = 2*diff(mnWav(2:end-1,:)) + (mnWav(4:end,:) - mnWav(1:end-3,:));
        mnWav([1,end-1,end],:) = 0;
        case 3
        mnWav(3:end-3,:) = 3*diff(mnWav(3:end-2,:)) + 2*(mnWav(5:end-1,:) - mnWav(2:end-4,:)) + (mnWav(6:end,:) - mnWav(1:end-5,:));
        mnWav([1,2,end-2,end-1,end],:) = 0;
        otherwise
        disperr_('ndiff_: nDiff_filt is valid for 1-3');
    end
    % catch
    %     if isGpu_(mnWav)
    %         mnWav = ndiff_(jrclust.utils.tryGather(mnWav), nDiff_filt);
    %     else
    %         disperr_('ndiff');
    %     end
    % end
end %func
