%--------------------------------------------------------------------------
function samplesPreview = loadFilePreview(fidBinary, P)
    % Peek ahead in the file given by fidBinary and rewind
    if P.nPaddingSamples > 0
        [samplesPreview, ~, sampleDims] = load_file_(fidBinary, P.nPaddingSamples, P);
        frewind_(fidBinary, sampleDims, P.dataType);
    else
        samplesPreview = [];
    end
end % function
