classdef ndiRecording < jrclust.interfaces.RawRecording
    %NDIRECORDING Model of an element in NDI
    %% NDI-SPECIFIC PROPERTIES

    properties (GetAccess=public,SetAccess=private)
        S; % ndi.session.dir that holds the experimental session data
        E; % the ndi.element object that holds the data
        epoch_id; % the epoch ids that are included
        nsamples; % the number of samples per epoch
        t0_t1; % start and end times, in local device time, of each epoch
    end;

    %% LIFECYCLE    
    methods
        function obj = ndiRecording(epochname, hCfg)
            %NDIRECORDING Construct an instance of this class
            % set object data type
            obj = obj@jrclust.interfaces.RawRecording(epochname, hCfg);
            % file not found error is not an error for ndi
            obj.errMsg = '';
            obj.isError = 0;
            obj.S = ndi.session.dir(hCfg.ndiPath);
            E = getelements(obj.S,'element.name',hCfg.ndiElementName,'element.reference',hCfg.ndiElementReference);
            if iscell(E) & numel(E)==1,
                obj.E = E{1};
            else,
                obj.isError = 1;
                obj.errMsg = 'No such element found.';
                return;
            end;

            try
                obj.hCfg = hCfg;
            catch ME
                obj.errMsg = ME.message;
                obj.isError = 1;
                return;
            end

            et = epochtable(obj.E);
            
            [parentdir,epoch_id] = fileparts(epochname);
            match = find(strcmp(epoch_id,{et.epoch_id}));
            if isempty(match),
                 obj.isError = 1;
                 obj.errMsg = ['No such epoch ' epoch_id '.'];
                 return;
            end;
            obj.epoch_id = et(match).epoch_id;
            % find clock that is dev local time
            foundMatch = 0;
            for i=1:numel(et(match).epoch_clock),
                 if strcmp(et(match).epoch_clock{i}.type,'dev_local_time'),
                     foundMatch = i;
                     break;
                 end;
            end;
            if foundMatch,
                obj.t0_t1 = et(match).t0_t1{foundMatch};
            else,
                error(['No clock type ''dev_local_time''.']);
            end;
            obj.nsamples = 1+diff(times2samples(obj.E,epoch_id,obj.t0_t1));
            obj.dshape = [hCfg.nChans, obj.nsamples];
        end % ndiRecording()

        function openRaw(obj)
            %OPENRAW Open the raw recording file for reading - does nothing for ndiRecording
            obj.rawIsOpen = 1;
        end
        function closeRaw(obj)
            %CLOSERAW Close the raw recording file for reading - does nothing for ndiRecording
            obj.rawIsOpen = 0;
        end

        function roi = readRawROI(obj, rows, cols)
            %READRAWROI Get a region of interest by rows/cols from the raw file
            t0t1 = samples2times(obj.E, obj.epoch_id, cols([1 end]));
            roi = readtimeseries(obj.E, obj.epoch_id, t0t1(1), t0t1(2));
            roi = roi'; % switch to column-based samples
            roi = single(roi(rows,:)); % if only a subset requested, return only the subset
        end % readRawROI()


    end % methods

    %% GETTERS/SETTERS
    methods
        function set.S(obj, val)
            if ~isa(val,'ndi.session.dir'),
                 error(['S must be of type ndi.session.dir']);
            end;
            obj.S = val;
        end
        function set.E(obj,val)
            if ~isa(val,'ndi.element') & ~isa(val,'ndi.time.timeseries'),
                 error(['E must be an ndi.element and an ndi.time.timeseries']);
            end;
            obj.E = val;
        end;
	
    end
end

