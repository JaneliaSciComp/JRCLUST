% test functions

%% No parameters
S_out = jrc3a('test', 'help_', {}, 0)

S_out = jrc3a('test', 'doc_', {}, 0)

S_out = jrc3a('test', 'disperr_', {}, 0)

jrc3a('test', 'setpath_', {}, 0)
clear jrc3a
jrc3a('test', 'setpath_', {}, 0)
jrc3a('test', 'setpath_', {}, 0)

jrc3a('test', 'install_', {}, 0)

jrc3a('test', 'disp_dependencies_', {}, 0)

jrc('test', 'check_requirements_', {}, 0)

csFiles = {'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_behav_g0_t1*.bin', 'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_behav_g0_t2*.bin'};
S_out = jrc3('test', 'filter_files_', {csFiles}, 2)

csFiles = {'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_behav_g0_t1*.bin', 'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_behav_g0_t0.nidq.bin'};
S_out = jrc3('test', 'filter_files_', {csFiles}, 2)

csFiles = {'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_behav_g0_t1*.bin', 'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_behav_g0_t0.nidq.bin', 'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_behav_g0_t2*.bin'};
S_out = jrc3('test', 'filter_files_', {csFiles}, 2)

S_out = jrc3('test', 'dir_file_', {'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_behav_g0_t0.nidq.bin', 1}, 2)

S_out = jrc3('test', 'list_files_', {{'sample*.bin'}}, 1)
S_out = jrc3('test', 'list_files_', {{'sample.bin', 'sample3.bin'}}, 1)
S_out = jrc3('test', 'list_files_', {{'sample*.bin', 'sample3.bin'}}, 1)


S_out = jrc3('test', 'subsFileExt_', {'c:\test.txt', '.prm'}, 1)
S_out = jrc3('test', 'subsFileExt_', {'c:\test.txt', '.prm', '.exe'}, 2)

S_out = jrc3('test', 'struct_merge_', {struct('a',1), struct('b',2,'c',3)}, 1)
S_out = jrc3('test', 'struct_merge_', {struct('b',1), struct('b',2,'c',3)}, 1)

S_out = jrc3('test', 'file_info_', {'sample.bin'}, 1);
S_out = jrc3('test', 'file_info_', {'sample.bi'}, 1);

S_out = jrc3('test', 'clear_', {}, 0);
S_out = jrc3('test', 'clear_', {'sample_sample.prm'}, 0);

S_out = jrc3('test', 'clear_', {}, 0);

S_out = jrc3('test', 'delete_empty_files_', {}, 0);
S_out = jrc3('test', 'find_empty_files_', {}, 1);

copyfile jrc3.m jrc3b.m f
S_out = jrc3('test', 'delete_files_', {'jrc3b.m', 1}, 0);
S_out = jrc3('test', 'delete_files_', {'jrc3c.m', 1}, 0);

S_out = jrc3a('test', 'doc_', {}, 0);

%% Test makeprm command
S_out = jrc3('test', 'makeprm_', {'sample.bin', 'sample.prb'}, 1);
S_out = jrc3('test', 'read_meta_file_', {'sample.meta'}, 1);
S_out = jrc3('test', 'read_meta_file_', {'sample.met'}, 1);
S_out = jrc3('test', 'read_meta_file_', {'sample3.meta'}, 1);
S_out = jrc3('test', 'read_whisper_meta_', {'sample3.meta'}, 1);
S_out = jrc3('test', 'text2struct_', {'sample3.meta'}, 1);

S_out = jrc3('test', 'get0_', {''}, 1);

S_out = jrc3('test', 'assignWorkspace_', {''}, 0);

S_out = jrc3('test', 'bytesPerSample_', {'int16'}, 1);
S_out = jrc3('test', 'bytesPerSample_', {'double'}, 1);
S_out = jrc3('test', 'bytesPerSample_', {'single'}, 1);
S_out = jrc3('test', 'bytesPerSample_', {'albert'}, 1);

S_out = jrc3a('test', 'makeprm_', {'sample.bin', 'sample.prb'}, 1);
S_out = jrc3('test', 'read_meta_file_', {'sample.meta'}, 1);
S_out = jrc3('test', 'read_meta_file_', {'sample.met'}, 1);
S_out = jrc3('test', 'read_meta_file_', {'sample3.meta'}, 1);
S_out = jrc3('test', 'read_whisper_meta_', {'sample3.meta'}, 1);
S_out = jrc3('test', 'text2struct_', {'sample3.meta'}, 1);
S_out = jrc3('test', 'get0_', {''}, 1);
S_out = jrc3('test', 'assignWorkspace_', {''}, 1);

S_out = jrc3('test', 'first_string_', {'e = mc^2; % hello'}, 1);
S_out = jrc3('test', 'first_string_', {' e = mc^2; % hello'}, 1);

S_out = jrc3('test', 'file2cellstr_', {'default.prm'}, 1);
S_out = jrc3('test', 'getCommentExpr_', {'e = mc^2; % hello'}, 1);

a = 1; S_out = jrc3('test', 'field2str_', {a}, 1);
a = [2 3 4]; S_out = jrc3('test', 'field2str_', {a}, 1);
a = [2 3.1 4]; S_out = jrc3('test', 'field2str_', {a}, 1);
a = {'a', 'b', 'c'}; S_out = jrc3('test', 'field2str_', {a}, 1);
a = {'a', 'b', 3}; S_out = jrc3('test', 'field2str_', {a}, 1);
a = struct('a',1,'b',2); ; S_out = jrc3('test', 'field2str_', {a}, 1);

S_out = jrc3('test', 'cellstr2file_', {'test.txt', {'a', 'b=c;', 'd=''efg'''}}, 0); edit('test.txt');

a=1; S_out = jrc3('test', 'set0_', {a}, 1);


%% 
S_out = jrc3a('test', 'read_cfg_', {}, 1)

S_out = jrc3a('test', 'makeStruct_', {}, 1)

S_out = jrc3a('test', 'matchFileExt_', {'trialstamp.csv', '.csv'}, 1)
S_out = jrc3a('test', 'matchFileExt_', {'trialstamp.cv', '.csv'}, 1)
S_out = jrc3a('test', 'matchFileExt_', {{'trialstamp.csv', 'test.cv', 'test.csv'}, '.csv'}, 1)

S_out = jrc3a('test', 'loadTrial_', {'trialstamp.csv'}, 1)

S_out = jrc3a('test', 'toString_', {112, 113}, 1)
S_out = jrc3a('test', 'toString_', {{'a','b','c'}}, 1)

S_out = jrc3a('test', 'squeeze_', ones(1,10,3))
S_out = jrc3a('test', 'squeeze_', ones(1,10,3), 1)
S_out = jrc3a('test', 'squeeze_', {ones(1,10,3)}, 1)

S_out = jrc3('test', 'signlog_', [-exp(1),exp(2),exp(3)])

S_out = jrc3('test', 'signsqrt_', [-2^2,3^2,4^2])


%% 8/4/17 JJJ: match file ext
S_out = jrc3('test', 'isTextFile_', {'silico_drift.Batch'}, 1)
S_out = jrc3('test', 'isTextFile_', {'silico_drift.txt'}, 1)
S_out = jrc3('test', 'isTextFile_', {'silico_drift.exe'}, 1)
S_out = jrc3('test', 'isTextFile_', {'silico_drift.exe', '.exe'}, 1);
S_out = jrc3('test', 'isTextFile_', {'silico_drift.exe', '.txt'}, 1);
S_out = jrc3('test', 'isTextFile_', {'silico_drift.exe', {'.txt', '.exe'}}, 1);

S_out = jrc3('test', 'filter_files_', {'silico_drift.batch'}, 2)

S_out = jrc3('test', 'batch_verify_', {'silico_drift.batch'}, 0)

S_out = jrc3('test', 'load_batch_', {'silico_drift.batch'}, 1);

S_out = jrc3('test', 'filter_files_', {'silico_drift.batch'}, 2)

S1 = struct('a', 1, 'b', 'hello');
S2 = struct('b', 2, 'c', 'cello');
S_out = jrc3('test', 'struct_merge_', {S1, S2}, 1);
S_out = jrc3('test', 'struct_merge_', {S1, S2, {'b'}}, 1);
S_out = jrc3('test', 'struct_merge_', {S1, S2, 'b'}, 1);

S_out = jrc3('test', 'read_cfg_', {}, 1);


%% 8/5/17 JJJ: subsample file
a = 1:10; S_out = jrc3('test', 'subsample_', {a, 3}, 1);
a = {'a', 'b','c','d'}; S_out = jrc3('test', 'subsample_', {a, 2}, 1);
a = {'a'}; S_out = jrc3('test', 'subsample_', {a, 30}, 1);

S_out = jrc3('test', 'sample_skip_', {[1 10], 100, 3}, 3)

S_out = jrc3('test', 'load_preview_', {'sample_sample.prm'}, 2, 0);
S_out1 = jrc3('test', 'load_preview_', {'sample_sample.prm'}, 1, 0);

a = rand(10,3);
S_out1 = jrc3('test', 'meanSubt_', {a}, 1, 0);
S_out2 = jrc3('test', 'meanSubt_', {a,1}, 1, 0);
S_out3 = jrc3('test', 'meanSubt_', {a,2}, 1, 0);
mean(S_out1.out1,1)
mean(S_out2.out1,1)
mean(S_out3.out1,1)
mean(S_out3.out1,2)
S_out4 = jrc3('test', 'meanSubt_', {a,1,@median}, 1, 0); median(S_out4.out1,1)

S_out = jrc3('test', 'isa_', {gpuArray(int16(rand(3,3))), 'int16'}, 1);
S_out = jrc3('test', 'isa_', {int16(rand(3,3)), 'int16'}, 1);
S_out = jrc3('test', 'isa_', {int16(rand(3,3)), 'single'}, 1);
S_out = jrc3('test', 'isa_', {single(rand(3,3)), 'single'}, 1);

hFig = figure('Tag', 'test');
S_out = jrc3('test', 'get_fig_cache_', {'test'}, 1);
S_out = jrc3('test', 'get_fig_cache_', {'test'}, 1);
S_out = jrc3('test', 'get_fig_cache_', {'tet'}, 1);

S_out = jrc3('test', 'preview_', {'sample_sample.prm'}, 1);


%% 8/7/17 JJJ: load batch

S_out = jrc3('test', 'load_batch_', {'silico_static.batch'}, 1);
jrc dir 'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\*.bin' mylist.txt

S_out = jrc3('test', 'file_created_', {'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_tag2_g0_t21.nidq.bin'}, 1);
datestr(S_out.out1, 'yyyy-mm-dd HH:MM:SS.FFF')
csFiles = {'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_tag2_g0_t21.nidq.bin', 'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_tag2_g0_t22.nidq.bin'};
S_out = jrc3('test', 'file_created_', {csFiles}, 1);
datestr(S_out.out1, 'yyyy-mm-dd HH:MM:SS.FFF')

vcFile_dir = 'E:\MikeEconomo\anm369962 - WR29\2017-03-06\SpikeGL\anm369962_*.nidq.bin';
S_out = jrc3('test', 'dir_file_', {vcFile_dir}, 1);


%% 8/15/17 JJJ: preview menu
S_out = jrc3('test', 'expand_vr_', {[1 3 5 7], 3}, 1);
S_out = jrc3('test', 'expand_vr_', {[1 3 5 7], 3, [11, 1]}, 1);
S_out = jrc3('test', 'expand_vr_', {[1 3 5 7], 3, [13,1]}, 1);


%% 9/6/17 JJJ: detrend funciton

S_out = jrc3('test', 'detrend_zscore_', {1:10, 1:10}, 2);
S_out = jrc3('test', 'detrend_zscore_', {[1:10], [1:9,0], 1:9}, 2);

S_out = jrc3('test', 'find_topn_', {[10 9 4 3 8 1 2], 3}, 1);
S_out = jrc3('test', 'find_topn_', {[10 9 4 3 8 1 2], 3, [1 2 3 4 6 7]}, 1);

S_out = jrc3('test', 'cell2vec_', {{[1,2,3], [4,5], [], [6,7]}}, 1);
S_out = jrc3('test', 'cell2vec_', {{[1,2,3], [4,5]', [], [6,7]}}, 1);



%%
S_out = jrc3('test', 'setlim_', {1:5, [3 5]}, 1);
S_out = jrc3('test', 'setlim_', {rand(3,3), [.3 .5]}, 1);


%% 9/26/17 JJJ
S_in = struct('a',1,'b',2,'c',struct('a',1,'b',2));
S_out = jrc3('test', 'struct_get_', {S_in, 'a','b','c','d'}, 1);

% Testing copy file
S_out = jrc3('test', 'copyfile_', {'c:\test1\*', 'c:\test2\'}, 0);
S_out = jrc3('test', 'copyfile_', {'c:\test1\*', {'c:\test2\', 'c:\test3\'}}, 0);
S_out = jrc3('test', 'copyfile_', {{'c:\test1\file1.txt', 'c:\test1\file2.txt'}, 'c:\test2\'}, 0);

S_out = jrc3('test', 'subsDir_', {'c:\test1\test.prb', 'd:\github\'}, 1);
S_out = jrc3('test', 'subsDir_', {'c:\test1\test.prb', 'd:\github'}, 1);
S_out = jrc3('test', 'subsDir_', {'c:\test1\test.prb', './probe'}, 1);

S_out = jrc3('test', 'search_file_', {'sample.prb', './prb'}, 1);
S_out = jrc3('test', 'search_file_', {'sample.prb', {'.pra', './prb'}}, 1);
S_out = jrc3('test', 'search_file_', {'sample.prb', {'.pra', './prc'}}, 1);

S_out = jrc3('test', 'jrcpath_', {}, 1)
S_out = jrc3('test', 'jrcpath_', {'test.m'}, 1)

S_out = jrc3('test', 'exist_file_', {'sample.bin'}, 1)
S_out = jrc3('test', 'exist_file_', {'./jrclust_alpha/sample.bin'}, 1)

%% 9/27/17 JJJ
S_out = jrc3('test', 'get_filter_', {struct('nDiff_filt', 0)}, 1)
S_out = jrc3('test', 'get_filter_', {struct('nDiff_filt', 1)}, 1)
S_out = jrc3('test', 'get_filter_', {struct('vcFilter', 'bandpass')}, 1)

S_out = jrc3('test', 'wiki_', {}, 0)
S_out = jrc3('test', 'wiki_', {'Command-line-interface'}, 0)

%% 9/29/17 JJJ
S_out = jrc3('test', 'jrc_version_', {}, 0);
S_out = jrc3('test', 'jrc_version_', {}, 1);
S_out = jrc3('test', 'jrc_version_', {}, 2);
S_out = jrc3('test', 'jrc_version_', {'sample_sample.prm'}, 2);
S_out = jrc3('test', 'jrc_version_', {'sample_sample.prm'}, 1);

%% 10/9/17
S_out = jrc3('test', 'importTSF_', {'E:\CatalinMitelut\sep_3_2017\layout_0.tsf'}, 2);

%% 10/11
S_out = jrc3('test', 'load_spkraw_', {}, 1, 0);
S_out = jrc3('test', 'load_spkwav_', {}, 1, 0);

S_out = jrc3('test', 'fread_spkwav_', {}, 1, 0);

S_out = jrc3('test', 'fread_spkwav_', {1:10}, 1, 0);
global tnWav_spk
b = tnWav_spk(:,:,1:10) - S_out.out1;
max(b(:)) - min(b(:))

S_out = jrc3('test', 'fread_spkraw_', {1:10}, 1, 0);
global tnWav_raw
b = tnWav_raw(:,:,1:10) - S_out.out1;
max(b(:)) - min(b(:))

stats = [];
viSpk = [1:100:190918];
viSpk2 = [viSpk(1:2:end); fliplr(viSpk(2:2:end))];
viSpk2 = viSpk2(:);
for i=1:30
    tic
    S_out = jrc3('test', 'fread_spkraw_', {viSpk}, 1, 0); %
    t1=toc;
    tic
%     global tnWav_spk
%     a = tnWav_spk(:,:,viSpk);
    S_out = jrc3('test', 'fread_spkraw_', {viSpk2}, 1, 0); %
    t2=toc;
    stats(end+1)=t1/t2;
end
printstat(stats)

%% 10/12

S_out = jrc3('test', 'set_diag_', {rand(3,4), nan}, 1); %
S_out = jrc3('test', 'set_diag_', {rand(3,3), nan}, 1); %
S_out = jrc3('test', 'set_diag_', {[], nan}, 1); %
S_out = jrc3('test', 'set_diag_', {rand(3,3), 1:3}, 1); %
S_out = jrc3('test', 'set_diag_', {rand(3,4), 1:3}, 1); %
S_out = jrc3('test', 'set_diag_', {rand(3,4), 1:4}, 1); %
