function bootstrapNDI(obj, varargin)
    %BOOTSTRAPNDI Build a configuration PRM file for reading an ndi.element.timeseries object
    %
    % BOOTSTRAPNDI(OBJ, NDI_SESSION_OBJ, NDI_ELEMENT_TIMESERIES_OBJ ...])
    %
    % Generates a PRM file for the ndi.element.timeseries object specified that is part of
    % the ndi.session that is indicated. 
    %
    % The JRCLUST analysis files will be installed at the path to NDI_SESSION_DIR_OBJ in a folder
    % called '.JRCLUST' and another subfolder with the string name of NDI_ELEMENT_TIMESERIES_OBJ. 
    % 
    % During bootstrap, all epochs of NDI_ELEMENT_TIMESERIES_OBJ are added for extraction. They
    % can be edited down as needed in the editor.
    % 
    % Example:
    % 
    % 

    if numel(varargin)~=2,
       error(['Usage: jrc(''bootstrap'',''ndi'',ndi_session_dir_obj,ndi_element_timeseries_obj)']);
    end;

    % Step 1: check our arguments and set our variables
    %   S - the ndi.session.dir
    %   E - the ndi.element.timeseries
    %   epoch_ids - the epoch ids

    try,
       ndi_v = ndi.version();
    catch,
       error(['Could not find ndi. Need to install https://github.com/VH-Lab/NDI-matlab ']);
    end;

    % confirm inputs
    try,
       S = varargin{1}; 
       if ~isa(S,'ndi.session.dir'),
           error(['Input argument that follows NDI is not of type ndi.session.dir']);
       end;
       E = varargin{2};
       if ~isa(E,'ndi.element') | ~isa(E,'ndi.time.timeseries'),
           error(['Input that follows ndi.session is not of type ndi.element.timeseries.']);
       end;
       if E.session~=S,
           error(['The ndi.element.timeseries must be part of the ndi.session provided.']);
       end;
       et = E.epochtable();
       epoch_ids = {et.epoch_id};
    catch ME,
       disp(['Usage: jrc(''bootstrap'',''ndi'',ndi_session_dir_obj,ndi_element_timeseries_obj)']);
       error(ME.message);
    end;

    Estring = E.elementstring();
    Estring(find(Estring==' ')) = '_';
    output_dir = [S.path() filesep '.JRCLUST' filesep Estring];
    output_file = [output_dir filesep 'jrclust.prm'];
    if ~exist(output_dir,'dir'),
        try,
           mkdir(output_dir);
        end;
    end;

    cfgData = struct();

    [data,t,timeref] = E.readtimeseries(1,0,0); % read 1 sample
    cfgData.tallSkinny = 1; % we will transpose later so this is met
    cfgData.outputDir = output_dir;

    cfgData.headerOffset = 0; % does not make sense for NDI
    cfgData.nChans = size(data,2);
    cfgData.sampleRate = E.samplerate(1); % get sample rate of first epoch; assume it's the same
    cfgData.rawRecordings = epoch_ids;
    
    cfgData.dataTypeRaw = 'int16'; % meaningless
    cfgData.dataTypeExtracted = 'single'; %  we will convert to single resolution
    cfgData.bitScaling = 1;
    cfgData.recordingFormat = 'ndi';
    cfgData.ndiPath = S.path();
    cfgData.ndiElementName = E.name;
    cfgData.ndiElementReference = E.reference;

    q = ndi.query('','depends_on',E.id(),'') & ndi.query('','isa','probe file','');
    docs = S.database_search(q);
    if ~isempty(docs),
         error(['I do not know how to do this yet.']);
         % set probe parameters from doc
    else,
            % set some trivial defaults
            cfgData.siteMap = 1:cfgData.nChans;
            cfgData.siteLoc = [zeros(cfgData.nChans, 1), 50*(0:cfgData.nChans-1)'];
            cfgData.shankMap = ones(cfgData.nChans, 1);
    end

    % construct the Config object from specified data
    hCfg_ = jrclust.Config(cfgData);

    configFile = fullfile(hCfg_.outputDir, ['jrclust', '.prm']);
    if ~exist(configFile,'file'),
        fclose(fopen(configFile,'w')); % touch the file and close it
    end;
    hCfg_.setConfigFile(configFile, 0);
    hCfg_.save('', 1);

    obj.hCfg = hCfg_;
    obj.hCfg.edit();
end

