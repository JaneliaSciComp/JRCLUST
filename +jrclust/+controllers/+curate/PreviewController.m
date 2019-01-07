classdef PreviewController
    %PREVIEWCONTROLLER
    
    properties (SetAccess=private)
        fileBounds;
        hCfg;
    end
    
    %% LIFECYCLE
    methods
        function obj = PreviewController(hCfg)
            %PREVIEWCONTROLLER Construct an instance of this class
            obj.hCfg = hCfg;
            obj.fileBounds = containers.Map();
        end
    end

    %% UTILITY METHODS
    methods (Hidden)
        function loadPreview(obj)
            nLoadsMax = obj.hCfg.nLoads_max_preview;
            nSecsLoad = obj.hCfg.sec_per_load_preview;

            rawRecordings = jrclust.utils.subsample(obj.hCfg.rawRecordings, nLoadsMax);

            % load files
            nLoadsFile = floor(nLoadsMax / numel(rawRecordings));
            nSamplesLoad = ceil(nSecsLoad*obj.hCfg.sampleRate);

            obj.hCfg.useGPU = false;

            for iFile = 1:numel(rawRecordings)
                hRec = jrclust.models.recording.Recording(rawRecordings{iFile}, obj.hCfg);
                if hRec.nSamples < nSamplesLoad
                    nLoadsiFile = 1;
                    nSamplesiLoad = hRec.nSamples;
                else
                    nLoadsiFile = min(nLoadsFile, floor(hRec.nSamples / nSamplesLoad));
                    nSamplesiLoad = nSamplesLoad;
                end

                multiBounds = sample_skip_([1, nSamplesiLoad], hRec.nSamples, nLoadsiFile);
                obj.fileBounds(hRec.binpath) = multiBounds;
            end
        end
    end
end

