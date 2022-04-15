function ndi(hCfg, varargin)
    % NDI - export spike clusters to NDI
    % 
    % NDI(HCFG, ...)
    %
    % Given a JRCLUST/Config object, export sorted spike clusters to NDI.
    %
    % The NDI session and NDI element are specified in the configuration information:
    % hCfg.ndiPath is the path to the NDI session to be modified.
    % hCfg.ElementName is the element name, and hCfg.ElementReference is the element reference.
    % The new neurons will be named [hCfg.ElementName '_N'], where N is the cluster number.
    % 
    % An NDI document of type JRCLUST_parameters is created and added to the NDI session.
    % It depends on the NDI element being sorted, and has a name equal to the JRCLUST parameter file name
    % (usually just 'jrclust.prm' unless the user has changed it). Its binary component is the text of the
    % JRCLUST parameter file. It also has a field JRCLUST.checksum that has the MD5 checksum of the
    % the jrclust_res.mat file. This checksum is used to check if the clustering has changed since a previous
    % export.
    %
    % If the clusters have not changed since a previous writing, no action is taken unless
    % the 'force' parameter is set to 1. If the clusters have changed, then the user is asked to
    % confirm that he/she wants to delete the previous cluster document (and the neurons that depend on
    % it) before re-writing.
    % 
    % This function takes additional name/value pairs that modify its operation:
    %
    % Parameter (default)    | Description
    % ----------------------------------------------------------------
    % interactive (1)        | Ask before deleting existing clusters
    % interactviakey         | Interact via keyboard (otherwise, use dialog)
    %   (hCfg.batch)         |    
    % forceReplace (0)       | Re-write clusters even if they have not changed
    %                        |   (will delete any documents that depend on the neurons,
    %                        |   including analyses.)
    %

        % Step 1a: Establish default parameters:

    interactive = 1;
    interactiveviakey = hCfg.batchMode;
    forceReplace = 0;

        % Step 1b: Update parameters with any user requests
    jrclust.utils.assign(varargin{:});

        % Step 1c: Check that we have NDI
    try,
        ndi_v = ndi.version();
    catch,
        error(['Could not determine NDI version. Check that https://github.com/VH-Lab/NDI-matlab is installed.']);
    end;

       % Step 2: load NDI session and element and check for errors

    disp(['Opening NDI experiment at ' hCfg.ndiPath '...']);
 
    S = ndi.session.dir(hCfg.ndiPath);
    E = getelements(S,'element.name',hCfg.ndiElementName,'element.reference',hCfg.ndiElementReference);
    if iscell(E) & numel(E)==1,
        E = E{1};
    else,
        error(['No such element ' hCfg.ndiElementName ' with reference ' num2str(hCfg.ndiElementReference) ' found.']);
    end;

       % Step 3: check for existing cluster document; if it exists, ask if we want to delete it

    md5_value = jrclust.utils.DataHash.DataHash(hCfg.resFile,'file','MD5');

    [configPath, configFileName, configFileExt] = fileparts(hCfg.configFile);
    configFileName = [configFileName configFileExt];
   
    q = ndi.query('','depends_on','element_id',E.id()) & ...
       ndi.query('','isa','jrclust_clusters','') & ...
       ndi.query('ndi_document.name','exact_string', configFileName,'');
    doc = S.database_search(q);
    if isempty(doc), % no docs,
        % nothing to do
    elseif numel(doc)>1,
        % somehow we have more than one of these documents; this should not happen
        % we will gracefully fix it
        % 
        fixDuplicates = 1;
        if interactive,
            if interactiveviakey,
		disp(['We found more than 1 cluster document that matches this element and JRCLUST run (we found ' int2str(numel(doc)) ').']);
                r = input(['Shall we delete them all? (Cannot continue if not yes.) [Y/n] :']);
                if ~strcmpi(r,'Y') | ~strcmpi(r,'yes'),
                    fixDuplicates = 0;
                end;
            else,
                ButtonName = questdlg(['We found more than 1 cluster document that matches this element and JRCLUST run ' ...
                                       '(we found (' int2str(numel(doc)) '). ' ...
                                       'Shall we delete them all? (Cannot continue if not yes)'],'Delete extra cluster documents?',...
                                       'Yes','No/Cancel','Yes');
                if ~strcmp(ButtonName,'Yes'),
                    fixDuplicates = 0;
                end;
            end;
         end;
         if fixDuplicates,
            S.database_rm(doc);
            doc = {};
         else,
             error(['We cannot continue with erroneous duplicate documents.']);
         end;
    else, % we have 1 doc
        % see if it is up to date
        up2date = strcmp(doc{1}.document_properties.jrclust_clusters.res_mat_MD5_checksum,md5_value);
        if up2date,
            if ~forceReplace, 
               disp('we are already up to date, exiting...')
               return;
            else,
               S.database_rm(doc);
               doc = {};
            end;
        else,
            replaceIt = 1; % default will be to replace
            if ~forceReplace &interactive,
                if interactiveviakey,
                    r = input(['Update clusters to newer version held by JRCLUST? (Y/N): ']);
                    if ~strcmpi(r,'Y') | ~strcmpi(r,'yes'),
                        replaceIt = 0;
                    end;
                elseif ~forceReplace,
                    ButtonName = questdlg(['Update clusters to newer version held by JRCLUST?'],...
                                       'Update?', 'Yes','No/Cancel','Yes');
                    if ~strcmp(ButtonName,'Yes'),
                        replaceIt = 0;
                    end;
                end;
            end;
            if replaceIt, 
              S.database_rm(doc);
              doc = {};
            else, 
               return; % nothing to do
            end;
         end;
     end;
    
     % Step 4: If we are continuing, write the cluster document

     d = ndi.document('apps/JRCLUST/jrclust_clusters.json', ...
         'ndi_document.name', configFileName, ...
         'jrclust_clusters.res_mat_MD5_checksum', md5_value);
     d = d.set_dependency_value('element_id',E.id());
     S.database_add(d);

     % Step 5: Write each neuron; each should depend on the cluster document

     % Step 5a: load the spike definitions and make sure they are annotated
     c = load(hCfg.resFile);
     if ~isfield(c,'clusterNotes'),
         error(['Spikes have not yet been annotated. Annotate spikes first.']);
     end;

     % Step 5b: determine the mapping between samples and the NDI element's epoch

     sample_nums = [];
     epoch_ids = {};

     for i=1:numel(hCfg.rawRecordings),
         disp(['Examining epoch ' hCfg.rawRecordings{i} '...']);
         [epochpath, epochfile, epochext] = fileparts(hCfg.rawRecordings{i});
         epoch_ids{i} = [epochfile epochext];
         obj = jrclust.detect.ndiRecording(epoch_ids{i},hCfg);
         sample_nums(i) = obj.dshape(2);
         t0_t1s{i} = obj.t0_t1;
     end;

     % Step 5c: determine the clusters that are "good"

     clusters_to_output = [];
     for i=1:numel(c.clusterNotes),
        if ~strcmpi(c.clusterNotes{i},'noise'), % add other conditions here
            clusters_to_output(end+1) = i;
        end;
     end;

     sample_nums = [0 cumsum(sample_nums(:)') Inf];

     % Step 5d: write the clusters
     dependency = struct('name','jrclust_clusters_id','value',d.id());
         % assume all spikes are eligible to fire for all epochs

     matlab_ver = ver('MATLAB');
     matlab_version = matlab_ver.Version;
     app_struct = struct('name', 'JRCLUST', 'version', jrclust.utils.version, 'url', 'https://github.com/JaneliaSciComp/JRCLUST', ...
         'os', computer, 'os_version', '', 'interpreter', 'MATLAB', 'interpreter_version', matlab_version);

     for i=1:numel(clusters_to_output),
          element_neuron = ndi.neuron(S,[E.name '_' int2str(clusters_to_output(i))],...
		E.reference,'spikes',E,0,[],dependency);
          neuron_extracellular.number_of_samples_per_channel = size(c.meanWfGlobal,1);
          neuron_extracellular.number_of_channels = size(c.meanWfGlobal,2);
          neuron_extracellular.mean_waveform = squeeze(c.meanWfGlobal(:,:,clusters_to_output(i)));
          neuron_extracellular.waveform_sample_times= [hCfg.evtWindowSamp(1):hCfg.evtWindowSamp(2)] / hCfg.sampleRate;
          neuron_extracellular.cluster_index = clusters_to_output(i);
          if ischar(c.clusterNotes{clusters_to_output(i)})
              switch lower(c.clusterNotes{clusters_to_output(i)}),
                   case 'single', value = 1;
			disp(['Single unit cluster ' int2str(clusters_to_output(i)) ' (will be added): ' c.clusterNotes{clusters_to_output(i)} ]);
                   case 'multi', value = 4;
			disp(['Multi-unit cluster ' int2str(clusters_to_output(i)) ' (will be added): ' c.clusterNotes{clusters_to_output(i)} ]);
                   otherwise,
			disp(['Unknown cluster ' int2str(clusters_to_output(i)) ' (will be skipped): ' c.clusterNotes{clusters_to_output(i)} ]);
                       value = -1; % unsure
              end;
          else,
              disp(['Unknown cluster ' int2str(clusters_to_output(i)) ' (will be skipped): (empty note field)']);
              value = -1;
          end;
          if value<0, continue; end; % skip the cell if it is not even a multi-unit
          neuron_extracellular.quality_number = value;
          neuron_extracellular.quality_label = c.clusterNotes{clusters_to_output(i)};
          neuron_doc = ndi.document('neuron/neuron_extracellular.json','app',app_struct,'neuron_extracellular',neuron_extracellular);
          neuron_doc = neuron_doc.set_dependency_value('element_id',element_neuron.id());
          S.database_add(neuron_doc);
          spike_indexes = c.spikeTimes(c.spikesByCluster{clusters_to_output(i)});
          for j=1:numel(epoch_ids),
		if ~strcmp(epoch_ids{j},'t00002'),
                local_sample_indexes = find((spike_indexes > sample_nums(j)) & (spike_indexes <= sample_nums(j+1)));
                local_sample = spike_indexes(local_sample_indexes) - sample_nums(j); % convert to local samples
                spike_times_in_epoch = E.samples2times(epoch_ids{j},double(local_sample));
		element_neuron.addepoch(epoch_ids{j},ndi.time.clocktype('dev_local_time'),...
			t0_t1s{j},spike_times_in_epoch(:),ones(size(spike_times_in_epoch(:))));
		end;
          end;
     end;

end

