classdef Recording < handle & dynamicprops
    %RECORDING Model of a single recording
    properties (Hidden, SetAccess=protected, SetObservable)
        hCfg;           % Config object
    end

    properties (SetAccess=private, SetObservable, Hidden, Transient)
        errMsg;         % error message (if there is an error in the file)
        isError;        % flag, whether or not there is an error
        rawIsOpen;      % flag, whether raw recording file is open
        filtIsOpen;     % flag, whether filtered file is open
    end

    properties (SetAccess=private, SetObservable, Transient)
        rawData;        % raw data in memmapped form
        filteredData;   % filtered data in memmapped form
        filteredFid;    % file handle for writing filtered data
    end

    properties (SetObservable, Dependent, Transient)
        nChans;         % number of channels (rows) in the recording
        nSamples;       % number of samples per channel (columns) in the recording
    end

    properties (SetAccess=private, SetObservable)
        rawPath;        % absolute path to binary file
        filtPath;       % absolute path to filtered recording, if there is one
        metaPath;       % absolute path to meta file, if there is one

        startTime;      % beginning of recording, in samples
        endTime;        % end of recording, in samples

        dataType;       % data type contained in file
        dshape;         % shape of data (rows x columns), in samples
        fSizeBytes;     % size of the file, in bytes
        headerOffset;   % number of bytes at the beginning of the file to skip
    end

    %% LIFECYCLE
    methods
        function obj = Recording(filename, hCfg) % dataType, nChans, headerOffset, 
            %RECORDING Construct an instance of this class
            % check filename exists
            obj.rawPath = jrclust.utils.absPath(filename); % returns empty if not found
            obj.isError = isempty(obj.rawPath);
            if obj.isError
                obj.errMsg = sprintf('file not found: %s', filename);
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

            % set object data shape
            d = dir(obj.rawPath);
            obj.fSizeBytes = d.bytes;

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

            try
                obj.hCfg = hCfg;
            catch ME
                obj.errMsg = ME.message;
                obj.isError = 1;
                return;
            end

            obj.rawIsOpen = 0;
            obj.filtIsOpen = 0;
        end
    end

    %% GETTERS/SETTERS
    methods
        % hCfg
        function set.hCfg(obj, hc)
            failMsg = 'hCfg must be an instance of jrclust.Config';
            assert(isempty(hc) || isa(hc, 'jrclust.Config'), failMsg);
            obj.hCfg = hc;
        end

        % rawData
        function rd = get.rawData(obj)
            if obj.rawIsOpen
                rd = obj.rawData.Data.Data;
            else
                rd = [];
            end
        end

        % filteredData
        function val = get.filteredData(obj)
            if obj.rawIsOpen
                val = obj.filteredData.Data.Data;
            else
                val = [];
            end
        end

        % nChans
        function nc = get.nChans(obj)
            nc = obj.dshape(1);
        end

        % nSamples
        function nc = get.nSamples(obj)
            nc = obj.dshape(2);
        end
    end
end

