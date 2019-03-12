classdef Recording < jrclust.interfaces.RawRecording
    %RECORDING Model of a single SpikeGLX recording
    %% SPIKEGLX-SPECIFIC PROPERTIES
    properties (SetAccess=protected, SetObservable, Transient)
        rawData;        % raw data in memmapped form
    end

    properties (SetAccess=protected, SetObservable)
        startTime;      % beginning of recording, in samples
        endTime;        % end of recording, in samples
    end

    properties (SetObservable)
        metaPath;       % absolute path to meta file, if there is one
    end

    %% LIFECYCLE
    methods
        function obj = Recording(filename, hCfg) % dataType, nChans, headerOffset, 
            %RECORDING Construct an instance of this class
            % check filename exists
            obj = obj@jrclust.interfaces.RawRecording(filename, hCfg);
            if obj.isError
                return;
            end

            % set a filtered path
            [~, ~, ext] = fileparts(obj.rawPath);
            obj.filtPath = jrclust.utils.subsExt(obj.rawPath, ['.filtered' ext]);
            obj.filteredFid = -1;

            % set object data type
            obj.dataType = hCfg.dataType;

            % set headerOffset
            obj.headerOffset = hCfg.headerOffset;

            d = dir(obj.rawPath);
            nSamples = (d.bytes - obj.headerOffset) / hCfg.nChans / jrclust.utils.typeBytes(obj.dataType);
            if ceil(nSamples) ~= nSamples % must be an integer or we don't have a complete recording
                obj.errMsg = 'Number of samples computed is not an integer. Check your sample rate or nChans?';
                obj.isError = 1;
                return;
            end
            obj.dshape = [hCfg.nChans nSamples];

            % load start and end times
            obj.metaPath = jrclust.utils.findMeta(filename);
            if ~isempty(obj.metaPath)
                md = jrclust.utils.metaToStruct(obj.metaPath);
                if isfield(md, 'firstSample')
                    obj.startTime = md.firstSample;
                    obj.endTime = md.firstSample + obj.dshape(2);
                end
            end

            obj.rawIsOpen = 0;
            obj.filtIsOpen = 0;
        end
    end

    %% GETTERS/SETTERS
    methods
        % rawData
        function rd = get.rawData(obj)
            if obj.rawIsOpen
                rd = obj.rawData.Data.Data;
            else
                rd = [];
            end
        end
    end
end

