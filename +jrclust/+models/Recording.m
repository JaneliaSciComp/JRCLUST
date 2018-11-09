classdef Recording < handle
    %RECORDING Model of a single recording

    properties (Access=private, SetObservable, Hidden, Transient)
        fid;
        isMemmap;
        isOpen;
    end

    properties (SetAccess=private, SetObservable, Transient)
        fileData;
    end
    
    properties (SetObservable)
        filename;
    end

    methods
        function obj = Recording(varargin)
            %RECORDING Construct an instance of this class
        end
    end
end

