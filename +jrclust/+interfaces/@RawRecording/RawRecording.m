classdef (Abstract) RawRecording < handle
    %RAWRECORDING Model of a single recording
    properties (Hidden, SetObservable)
        hCfg;           % Config object
    end

    properties (Hidden, SetAccess=protected, SetObservable, Transient)
        rawIsOpen;      % flag, whether raw recording file is open
        filtIsOpen;     % flag, whether filtered file is open
    end

    properties (SetAccess=protected, SetObservable, Transient)
        errMsg;         % error message (if there is an error in the file)
        isError;        % flag, whether or not there is an error

        filteredData;   % filtered data in memmapped form
        filteredFid;    % file handle for writing filtered data
    end

    properties (Dependent, SetObservable, Transient)
        nChans;         % number of channels (rows) in the recording
        nSamples;       % number of samples per channel (columns) in the recording
    end

    properties (SetAccess=protected, SetObservable)
        rawPath;        % absolute path to binary file
        filtPath;       % absolute path to filtered recording, if there is one

        dataType;       % data type contained in file
        dshape;         % shape of data (rows x columns), in samples
        fSizeBytes;     % size of the file, in bytes
        headerOffset;   % number of bytes at the beginning of the file to skip
    end

    %% LIFECYCLE
    methods
        function obj = RawRecording(filename, hCfg)
            %RAWRECORDING Construct an instance of this class
            % check filename exists
            obj.rawPath = jrclust.utils.absPath(filename); % returns empty if not found
            obj.isError = isempty(obj.rawPath);
            if obj.isError
                obj.errMsg = sprintf('file not found: %s', filename);
                return;
            end

            try
                obj.hCfg = hCfg;
            catch ME
                obj.errMsg = ME.message;
                obj.isError = 1;
                return;
            end

            % set object data shape
            d = dir(obj.rawPath);
            obj.fSizeBytes = d.bytes;
        end
    end

    %% ABSTRACT METHODS
    methods (Abstract)
        roi = readRawROI(obj, rows, cols);
    end

    %% GETTERS/SETTERS
    methods
        % filteredData
        function val = get.filteredData(obj)
            if obj.rawIsOpen
                val = obj.filteredData.Data.Data;
            else
                val = [];
            end
        end

        % hCfg
        function set.hCfg(obj, hc)
            failMsg = 'hCfg must be an instance of jrclust.Config';
            assert(isempty(hc) || isa(hc, 'jrclust.Config'), failMsg);
            obj.hCfg = hc;
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

