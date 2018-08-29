%--------------------------------------------------------------------------
function mnWav1 = sgfilt_(mnWav, nDiff_filt, useGPU)
    % works for a vector, matrix and tensor
    fInvert_filter = 0;
    if nargin<3, useGPU = isGpu_(mnWav); end
    n1 = size(mnWav,1);
    if n1==1, n1 = size(mnWav,2);  end
    if nDiff_filt==0, mnWav1 = mnWav; return; end
    if fInvert_filter
        [miB, miA] = sgfilt_init_(n1, nDiff_filt, useGPU);
    else
        [miA, miB] = sgfilt_init_(n1, nDiff_filt, useGPU);
    end

    if isvector(mnWav)
        mnWav1 = mnWav(miA(:,1)) - mnWav(miB(:,1));
        for i=2:nDiff_filt
            mnWav1 = mnWav1 + i * (mnWav(miA(:,i)) - mnWav(miB(:,i)));
        end
    elseif ismatrix(mnWav)
        mnWav1 = mnWav(miA(:,1),:) - mnWav(miB(:,1),:);
        for i=2:nDiff_filt
            mnWav1 = mnWav1 + i * (mnWav(miA(:,i),:) - mnWav(miB(:,i),:));
        end
    else
        mnWav1 = mnWav(miA(:,1),:,:) - mnWav(miB(:,1),:,:);
        for i=2:nDiff_filt
            mnWav1 = mnWav1 + i * (mnWav(miA(:,i),:,:) - mnWav(miB(:,i),:,:));
        end
    end
end % function
