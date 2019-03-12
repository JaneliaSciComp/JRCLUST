function hRec = newRecording(rawPath, hCfg)
    %NEWRECORDING Create a new RawRecording based on recordingFormat
    switch hCfg.recordingFormat
        case 'Intan'
            hRec = jrclust.detect.IntanRecording(rawPath, hCfg);

        otherwise % SpikeGLX .bin/.dat
            hRec = jrclust.detect.Recording(rawPath, hCfg);

    end
end

