%--------------------------------------------------------------------------
function [cvnlim_bin, viRange, viEdges] = sample_skip_(nlim_bin, nSamples_bin, nTime_traces)
    % return a limit that
    % nlim_bin=[81 90]; nSamples_bin=100; nTime_traces=5;
    % edges to set to nan
    % 2017/6/22 James Jun: Added nTime_traces multiview
    % 6/23 JJJ: edge samples to be set to nan (gap)

    if nTime_traces==1 || isempty(nTime_traces)
        cvnlim_bin = {nlim_bin};
        viRange = nlim_bin(1):nlim_bin(end);
        viEdges = [];
        return;
    end
    nSkip = floor(nSamples_bin / nTime_traces);
    cvnlim_bin = arrayfun(@(i)nlim_bin + (i-1)*nSkip, 1:nTime_traces, 'UniformOutput', 0);
    % modulus algebra wrap around
    for i=1:nTime_traces
        lim1 = mod(cvnlim_bin{i}-1, nSamples_bin)+1;
        if lim1(1) > lim1(2)
            lim1 = [1, diff(nlim_bin)+1];
        end
        cvnlim_bin{i} = lim1;
    end
    if nargout>=2
        viRange = jrclust.utils.neCell2mat(cellfun(@(x)x(1):x(2), cvnlim_bin, 'UniformOutput', 0));
    end
    if nargout>=3 %compute the number of samples
        viEdges = cumsum(cellfun(@(x)diff(x)+1, cvnlim_bin));
        viEdges = [1, viEdges(1:end-1)];
        %     viEdges = sort([viEdges, viEdges+1], 'ascend'); %two sample gaps
    end
end %func
