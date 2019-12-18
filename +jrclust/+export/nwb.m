function nwbData = nwb(hClust, outPath)
    %NWBEXPORT exports jrclust data to NWB
    %   Assumes matnwb is in path and generateCore is called
    %   hClust is a regular hClust object with certain additions for config data
    %   outFileName is the destination file name produced by NWB.
    %   Should use the .nwb extension
    if exist('nwbRead', 'file') ~= 2
        warning('Please make sure you have MatNWB installed and on your path (https://github.com/NeurodataWithoutBorders/matnwb)');
        return;
    end

    assert(ischar(outPath), 'output filename should be a char array');

    if 2 == exist(outPath, 'file')
        [~, filename, ~] = fileparts(outPath);
        prompt = ['This file already exists.  '...
            'Would you like to overwrite or update the file?  '...
            'You will lose all data if you overwrite.  '...
            'Updating a NWB file is currently incomplete.'];
        title = 'Overwrite Prompt';
        selection = questdlg(prompt, title, 'Overwrite', 'Update', 'Cancel', 'Overwrite');

        switch selection
            case 'Overwrite'
                warning(toMsgID('File'),...
                    'Deleting `%s` and starting anew.', filename);
                delete(outPath);
            case 'Update'
                error(toMsgID('Internal', 'Todo'),...
                    'This feature is a work in progress.  Please wait warmly.');
            case {'Cancel' ''} % including X-ing out
                warning(toMsgID('Abort'), 'Aborting export.');
                return;
            otherwise
                error(toMsgID('Internal', 'InvalidSelection'),...
                    'Unknown questdlg selection %s', selection);
        end
    end

    hCfg = hClust.hCfg;

    %% Metadata
    log_progress('Filling out Metadata');
    nwbData = NwbFile;
    nwbData.identifier = hCfg.sessionName;
    nwbData.session_start_time = datetime(hClust.detectedOn, 'ConvertFrom', 'datenum');

    param_fid = fopen([hCfg.sessionName '-full.prm']);
    param_data = fread(param_fid, '*char');
    fclose(param_fid);
    nwbData.general_data_collection = param_data .';

    % TODO invalid data, replace with file create date if it exists.
    nwbData.timestamps_reference_time = datetime('now');

    nwbData.session_description = 'JRCLust Spike Data';
    nwbData.general_source_script_file_name = mfilename('fullpath');

    log_done();
    %% Devices
    log_progress('Filling out electrode hardware data');
    prompt = 'Please provide a probe name that was used to acquire this data.';
    probe_name = prompt_required_input(prompt);
    nwbData.general_devices.set(probe_name, types.core.Device());
    probe_path = ['/general/devices/' probe_name];
    %% Electrode Groups
    shank_ids = unique(hCfg.shankMap);
    num_shanks = length(shank_ids);
    groups_path_stub = '/general/extracellular_ephys/';
    Shank2Name_Map = containers.Map('KeyType', 'double', 'ValueType', 'char');
    for i = 1:num_shanks
        shank_number = shank_ids(i);
        name = sprintf('Shank%02d', shank_number);
        description = sprintf('shank %02d of %02d', shank_number, num_shanks);
        prompt = sprintf('What is the location of shank %02d?', shank_number);
        location = prompt_required_input(prompt);
        Electrode_Group = types.core.ElectrodeGroup(...
            'description', description,...
            'location', location,...
            'device', types.untyped.SoftLink(probe_path));
        nwbData.general_extracellular_ephys.set(name, Electrode_Group);
        Shank2Name_Map(shank_number) = name;
    end
    %% Electrode Sites
    %   Note: sites is a subset of channels.  They would be channels of interest.
    %   There is no location data for non-site channels, so we assume that all sites
    %   are valid electrodes
    num_sites = hCfg.nSites;
    vector_shape = [num_sites 1];
    site_locations = hCfg.siteLoc;
    invalid_str_column = {repmat({'N/A'}, vector_shape)}; % extra cell layer due to struct behavior
    zero_column = zeros(vector_shape);
    electrode_names = cell(vector_shape);
    electrode_references = types.untyped.ObjectView.empty;
    for i = 1:num_sites
        shank_name = Shank2Name_Map(hCfg.shankMap(i));
        electrode_names{i} = shank_name;
        path = [groups_path_stub shank_name];
        electrode_references(i) = types.untyped.ObjectView(path);
    end

    Electrode_Data = struct(...
        'x', site_locations(:, 1),...
        'y', site_locations(:, 1),...
        'z', zero_column,...
        'imp', zero_column,...
        'location', invalid_str_column,...
        'filtering', invalid_str_column,...
        'group', electrode_references,...
        'group_name', {electrode_names},...
        'channel_index', hCfg.siteMap);
    Descriptions = struct(...
        'channel_index', 'index to valid channels');
    nwbData.general_extracellular_ephys_electrodes =...
        to_electrode_table(Electrode_Data, Descriptions);
    electrodes_path = '/general/extracellular_ephys/electrodes';

    log_done();
    %% Recording
    % log_progress('Converting to HDF5 files and Linking (This will take a while)');
    % 
    % % calculate channel -> site indices
    % site_indices = zeros(Config.nChans, 1);
    % for i = 1:Config.nChans
    %     index = find(Electrode_Data.channel_index == i, 1);
    %     if ~isempty(index)
    %         site_indices(i) = index;
    %     end
    % end
    % Electrodes_Ref = to_table_region(site_indices,...
    %     'Sites associated with these channels (0 indicates invalid channel)',...
    %     electrodes_path);
    % 
    % % can be multiple raw files.  Assumed that they're in order.
    % raw_filenames = Config.rawRecordings;
    % start_sample = 0;
    % for i = 1:length(raw_filenames)
    %     raw = raw_filenames{i};
    %     File_Meta = dir(raw);
    %     sample_size  = Config.bytesPerSample * Config.nChans;
    %     num_samples = File_Meta.bytes / sample_size;
    %     data_shape = [Config.nChans num_samples];
    %     Raw_Series = types.core.ElectricalSeries(...
    %         'data', bin_to_hdf5(raw, data_shape, Config.dataType),...
    %         'starting_time', start_sample,...
    %         'starting_time_rate', Config.sampleRate,...
    %         'electrodes', Electrodes_Ref);
    %     nwb.acquisition.set(sprintf('RawRecording%02d', i), Raw_Series);
    %     start_sample = start_sample + num_samples;
    % end
    % 
    % log_done();
    %% Units
    log_progress('Organising spike units');

    num_clusters = length(hClust.spikesByCluster);
    cluster_shape = [num_clusters 1];
    spike_time_by_cluster = cell(cluster_shape);
    for i = 1:num_clusters
        spike_i = hClust.spikesByCluster{i};
        spike_time_by_cluster{i} = hClust.spikeTimes(spike_i);
    end

    means = squeeze(hClust.meanWfLocalRaw(:, 1, :));
    waveform_size = size(means, 1);
    num_waveforms = size(means, 2);
    means_per_cluster = mat2cell(means,...
        waveform_size,...
        ones(num_waveforms, 1));

    % determine group ref based on sites
    sites_per_cluster = hClust.clusterSites;
    electrodes_per_cluster = types.untyped.ObjectView.empty;
    Electrodes_Table = nwbData.general_extracellular_ephys_electrodes;
    groups_per_site = Electrodes_Table.vectordata.get('group').data;
    for i = length(sites_per_cluster)
        site_i = sites_per_cluster(i);
        electrodes_per_cluster(end+1) = groups_per_site(site_i);
    end

    Electrodes_Ref = struct(...
        'data', {num2cell(hClust.clusterSites)},...
        'table_path', electrodes_path);
    Cluster_Data = struct(...
        'spike_times', {spike_time_by_cluster},...
        'electrode_group', electrodes_per_cluster,...
        'electrodes', Electrodes_Ref,...
        'notes', {hClust.clusterNotes},...
        'waveform_mean', {means_per_cluster},...
        'peak_to_peak', hClust.unitVpp);
    if isprop(hClust, 'unitSNR')
        Cluster_Data.signal_to_noise = hClust.unitSNR;
    end

    Cluster_Descriptions = struct(...
        'spike_times', 'Spike times per cluster',...
        'electrode_group', 'shank used',...
        'electrodes', 'electrode site indices',...
        'notes', 'cluster-specific notes',...
        'waveform_mean', 'the spike waveform mean for each spike unit',...
        'signal_to_noise', 'signal to noise ratio for each unit.',...
        'peak_to_peak', 'Voltage peak to peak for this unit');
    nwbData.units = types.core.Units('description', 'Spikes organized by Cluster');
    nwbData.units = fill_dynamic_columns(nwbData.units, Cluster_Data, Cluster_Descriptions, '/units');

    log_done();
    %% WRITE
    log_progress('Writing NWB file');
    nwbExport(nwbData, outPath);
    log_done();
end

%% TO_INDEXED_VECTOR converts cell array data, description, and index path into a Vector[Index/Data] pair.
%   Return value is a struct containing 'data' -> VectorData and 'index' -> VectorIndex
function Vector = to_indexed_vector(data, description, index_path)
assert(iscell(data) && ~iscellstr(data), toMsgID('IndexedVector', 'InvalidArg'),...
    'To make a regular vector data object, use to_vector instead.');

indices = cumsum(cellfun('length', data));
ref = types.untyped.ObjectView(index_path);

is_vertical = all(cellfun('size', data, 1) > 1);
if is_vertical
    full_data = vertcat(data{:});
else
    full_data = [data{:}];
end
Vector = struct(...
    'data', to_vector(full_data, description),...
    'index', types.core.VectorIndex('data', indices, 'target', ref));
end

%% TO_VECTOR converts data and description to valid VectorData or Vector[Index/Data] pair.
function Vector_Data = to_vector(data, description)
assert(~iscell(data) || iscellstr(data), toMsgID('DataVector', 'InvalidArg'),...
    'To make an indexed vector data pair, use to_indexed_vector instead');
Vector_Data = types.core.VectorData('data', data, 'description', description);
end

%% TO_TABLE_REGION converts indices and h5_path to a DynamicTableRegion object
function Table_Region = to_table_region(indices, description, h5_path)
Table_Region = types.core.DynamicTableRegion(...
    'data', indices,...
    'description', description,...
    'table', types.untyped.ObjectView(h5_path));
end

%% TO_ELECTRODE_TABLE converts a struct of column names and data with column data.
function Electrode_Table = to_electrode_table(Column_Data, varargin)
% ripped from NWB format since electrode table isn't an actual type
Descriptions = struct(...
    'x', 'the x coordinate of the channel location',...
    'y', 'the y coordinate of the channel location',...
    'z', 'the z coordinate of the channel location',...
    'imp', 'the impedance of the channel',...
    'location', 'the location of channel within the subject e.g. brain region',...
    'filtering', 'description of hardware filtering',...
    'group', 'a reference to the ElectrodeGroup this electrode is a part of',...
    'group_name', 'the name of the ElectrodeGroup this electrode is a part of');
if ~isempty(varargin)
    Extra_Descriptions = varargin{1};
    extra_fields = fieldnames(Extra_Descriptions);
    for i = 1:length(extra_fields)
        field = extra_fields{i};
        Descriptions.(field) = Extra_Descriptions.(field);
    end
end

Electrode_Table = types.core.DynamicTable('description', 'relevent electrodes in experiment');
Electrode_Table = fill_dynamic_columns(Electrode_Table, Column_Data, Descriptions);
end

%% FILL_DYNAMIC_COLUMNS given a DynamicTable or a subclass of DynamicTable, fills available column data
function Table = fill_dynamic_columns(Table, Column_Data, Column_Description, varargin)
error_stub = 'FillColumns';
if isempty(varargin)
    ref_path = ''; 
else
    ref_path = varargin{1};
end
assert(ischar(ref_path), toMsgID(error_stub, 'InvalidArgs'),...
    'reference path must be a valid HDF5 path.');

column_names = fieldnames(Column_Data);
column_lengths = zeros(size(column_names));
for name_i = 1:length(column_names)
    name = column_names{name_i};
    data = Column_Data.(name);
    description = Column_Description.(name);
    
    is_table_region = isstruct(data);
    if is_table_region
        Region = data;
        data = Region.data;
    end
    
    [Data, Vector_Index] = compose_column_data(name, data, description, ref_path);
    
    if is_table_region
        Data = to_table_region(Data.data, description, Region.table_path);
    end
    
    assign_to_table(Table, name, Data, Vector_Index);
end
table_height = unique(column_lengths);
assert(isscalar(table_height), toMsgID(error_stub, 'InvalidHeight'),...
    'table columns must have the same length');

Table.colnames = column_names;
Table.id = types.core.ElementIdentifiers('data', 1:table_height);

    function [Data, Index] = compose_column_data(name, data, description, ref_path)
        is_indexed = iscell(data) && ~iscellstr(data);
        if is_indexed
            assert(~isempty(ref_path), toMsgID(error_stub, 'InvalidPath'),...
                'Indexed data detected in column but reference path was not provided.');
            full_ref_path = [ref_path '/' name];
            Vector = to_indexed_vector(data, description, full_ref_path);
            
            Data = Vector.data;
            Index = Vector.index;
        else
            Data = to_vector(data, description);
            Index = [];
        end
    end

    function assign_to_table(Table, name, data, varargin)
       if isprop(Table, name)
           Table.(name) = data;
       else
           Table.vectordata.set(name, data);
       end
       
       if ~isempty(varargin) && ~isempty(varargin{1})
           index = varargin{1};
           index_name = [name '_index'];
           if isprop(Table, index_name)
               Table.(index_name) = index;
           else
               Table.vectorindex.set(index_name, index);
           end
       end
    end
end

%% PROMPT_REQUIRED_INPUT creates a dialog asking the user to input data
%   The return value is the inputted string data.
function input = prompt_required_input(query)
input = inputdlg(query);
assert(~isempty(input),...
    toMsgID('Abort'),...
    'Required input was canceled.  Aborting export.');
input = input{1};
end

%% LOG_PROGRESS logs to command window.  Completed with LOG_DONE
function log_progress(log)
fprintf('%s...\n', log);
end

%% LOG_DONE returns Done after a LOG_PROGRESS call
function log_done()
fprintf('\bDone!\n');
end

%% TO_MSG_ID returns proper filterable msg_id used in error, warning, assert, etc.
function msgID = toMsgID(varargin)
    rootID = 'NwbExport';
    if isempty(varargin) || ~iscellstr(varargin)
        error(toMsgID('Internal', 'InvalidArg'),...
        'arguments must be non-empty cell strings.');
    end

    msgID = strjoin([{rootID} varargin], ':');
end

%% BIN_TO_HDF5 utility function to move binary data files to HDF5 equivalent.
function External_Link = bin_to_hdf5(filename, shape, datatype)
[path, name, ~] = fileparts(filename);
destination = [fullfile(path, name) '.h5'];

fcpl = H5P.create('H5P_FILE_CREATE');
fapl = H5P.create('H5P_FILE_ACCESS');
h5_fid = H5F.create(destination,'H5F_ACC_TRUNC', fcpl, fapl);
H5P.close(fcpl);
H5P.close(fapl);
fid = fopen(filename);

h5_type = mat2h5_datatype(datatype);
h5_shape = fliplr(shape);
space_id = H5S.create_simple(length(shape), h5_shape, h5_shape);
dataset_id = H5D.create(h5_fid, '/data', h5_type, space_id, 'H5P_DEFAULT');

h5_count = zeros(length(shape),1);
stripe_length = shape(1);
h5_count(end) = stripe_length;
for dim_i = 2:length(shape)
    factors = factor(shape(dim_i));
    largest_factor = max(factors);
    largest_selection = largest_factor^(sum(factors == largest_factor));
    stripe_length = stripe_length * largest_selection;
    h5_count(end - dim_i + 1) = largest_selection;
end

h5_start = zeros(length(shape), 1);
wait_bar = waitbar(0, 'Copying');
while ~feof(fid)
    data = fread(fid, stripe_length, ['*' datatype]);
    
    H5S.select_hyperslab(space_id, 'H5S_SELECT_SET', h5_start, [], h5_count, []);
    
    mem_space_id = H5S.create_simple(1, stripe_length, stripe_length);
    H5D.write(dataset_id, h5_type, mem_space_id, space_id, 'H5P_DEFAULT', data);
    H5S.close(mem_space_id);
    h5_start(1) = h5_start(1) + h5_count(1);
    waitbar(h5_start(1) / h5_shape(1), wait_bar,...
        sprintf('Copying %d/%d', h5_start(1), h5_shape(1)));
end

waitbar(1, 'Done');

H5S.close(space_id);
H5D.close(dataset_id);
H5F.close(h5_fid);
fclose(fid);
External_Link = types.untyped.ExternalLink(destination, '/data');

    function h5_datatype = mat2h5_datatype(mat_datatype)
        assert(startsWith(mat_datatype, {'int', 'uint'}),...
            'bin_to_h5 currently only supports integer data types');
        
        int_pattern = 'u?int(8|16|32|64)$';
        type_width = regexp(mat_datatype, int_pattern, 'tokens', 'once');
        assert(~isempty(type_width), 'invalid/unsupported datatype %s', mat_datatype);
        
        if mat_datatype(1) == 'u'
            unsigned_flag = 'U';
        else
            unsigned_flag = 'I';
        end
        
        h5_datatype = ['H5T_STD_' unsigned_flag type_width{1} 'LE'];
    end
end