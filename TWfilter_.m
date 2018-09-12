function TWfilter_(P, vcMode)
	% display manual sorting interface
	global fDebug_ui trFet_spk

	% Load info
	if ~is_sorted_(P)
		fprintf(2, 'File must to be sorted first (run "jrc spikesort %s")\n', P.vcFile_prm); 
		return; 
	end
	[S0, P] = load_cached_(P);
	if ~isfield(S0, 'mrPos_spk')
		S0.mrPos_spk = spk_pos_(S0, trFet_spk);
		set(0, 'UserData', S0);
	end
	
	fDebug_ui = 0;
	P.fGpu = 0; %do not use GPU for manual use
	set0_(fDebug_ui, P);

	S0 = set0_(P); %update the P structure
	S0.S_clu = S_clu_update_wav_(S0.S_clu, P);                
	set(0, 'UserData', S0);
    
    for i_clu = sort(find(S0.S_clu.viSite_clu>P.max_real_site),'descend')
    	S0.S_clu = delete_clu_(S0.S_clu, i_clu);
		fprintf('Deleting cluster %d\n', i_clu); 
    end
    set(0, 'UserData', S0);
	save0_(subsFileExt_(S0.P.vcFile_prm, '_jrc.mat'),1); % 1 will skip figure saving

end %func;

