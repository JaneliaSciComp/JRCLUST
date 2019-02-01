function samplesOut = sgFilter(samplesIn, nDiffOrder)
    %SGFILTER Savitzky-Golay filter
    % works for a vector, matrix and tensor
    fInvert_filter = 0;
    useGPU = isa(samplesIn, 'gpuArray');
    n1 = size(samplesIn,1);
    if n1 == 1
        n1 = size(samplesIn,2); 
    end

    if nDiffOrder == 0
        samplesOut = samplesIn;
        return;
    end

    [miA, miB] = sgfilt_init_(n1, nDiffOrder, useGPU);

    if isvector(samplesIn)
        samplesOut = samplesIn(miA(:,1)) - samplesIn(miB(:,1));
        for i = 2:nDiffOrder
            samplesOut = samplesOut + i * (samplesIn(miA(:,i)) - samplesIn(miB(:,i)));
        end
    elseif ismatrix(samplesIn)
        samplesOut = samplesIn(miA(:,1),:) - samplesIn(miB(:,1),:);
        for i = 2:nDiffOrder
            samplesOut = samplesOut + i * (samplesIn(miA(:,i),:) - samplesIn(miB(:,i),:));
        end
    else
        samplesOut = samplesIn(miA(:,1),:,:) - samplesIn(miB(:,1),:,:);
        for i = 2:nDiffOrder
            samplesOut = samplesOut + i * (samplesIn(miA(:,i),:,:) - samplesIn(miB(:,i),:,:));
        end
    end
end

%% LOCAL FUNCTIONS
function [miA, miB, viC] = sgfilt_init_(nData, nFilt, useGPU)
    persistent miA_ miB_ viC_ nData_ nFilt_
    if nargin<2, useGPU=0; end

    % Build filter coeff
    if isempty(nData_), nData_ = 0; end
    try a = size(miA_); catch, nData_ = 0; end
    if nData_ == nData && nFilt_ == nFilt
        [miA, miB, viC] = deal(miA_, miB_, viC_);
    else
        vi0 = jrclust.utils.tryGpuArray(int32(1):int32(nData), useGPU)';
        vi1 = int32(1):int32(nFilt);
        miA = min(max(bsxfun(@plus, vi0, vi1),1),nData);
        miB = min(max(bsxfun(@plus, vi0, -vi1),1),nData);
        viC = jrclust.utils.tryGpuArray(int32(-nFilt:nFilt), useGPU);
        [nData_, nFilt_, miA_, miB_, viC_] = deal(nData, nFilt, miA, miB, viC);
    end
end %func
