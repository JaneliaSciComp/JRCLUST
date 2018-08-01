%--------------------------------------------------------------------------
function nFailed = unit_test_(vcArg1, vcArg2, vcArg3)
    % 2017/2/24. James Jun. built-in unit test suite (request from Karel Svoboda)
    % run unit test
    %[Usage]
    % unit_test()
    %   run all
    % unit_test(iTest)
    %   run specific test again and show profile
    % unit_test('show')
    %   run specific test again and show profile
    % @TODO: test using multiple datasets and parameters.
    global fDebug_ui;

    if nargin<1, vcArg1 = ''; end
    if nargin<2, vcArg2 = ''; end
    if nargin<3, vcArg3 = ''; end

    cd(fileparts(mfilename('fullpath'))); % move to jrclust folder
    if ~exist_file_('sample.bin'), jrc3('download', 'sample'); end

    nFailed = 0;
    profile('clear'); %reset profile stats
    csCmd = {...
    'close all; clear all;', ... %start from blank
    'jrc3 clear sample_sample.prm', ...
    'jrc3 compile', ...
    'jrc3 probe sample.prb', ...
    'jrc3 makeprm sample_list.txt sample.prb', ...
    'jrc3 makeprm sample.bin sample.prb', ...
    'jrc3 probe sample_sample.prm', ...
    'jrc3 import-lfp sample_sample.prm', ...
    'jrc3 traces-lfp sample_sample.prm', ...
    'jrc3 traces sample_sample.prm', 'jrc3 traces', ...
    'jrc3 traces-test sample_sample.prm', ...
    'jrc3 preview-test sample_sample.prm', ...
    'jrc3 makeprm sample.bin sample.prb', ...
    'jrc3 detectsort sample_sample.prm', 'jrc3 clear sample_sample.prm', ...
    'jrc3 probe sample_sample_merge.prm', 'jrc3 detectsort sample_sample_merge.prm', 'jrc3 clear sample_sample_merge.prm', ... %multishank, multifile test
    'jrc3 detect sample_sample.prm', 'jrc3 sort sample_sample.prm', ...
    'jrc3 export-csv sample_sample.prm', ...
    'jrc3 export-quality sample_sample.prm', ...
    'jrc3 export-spkwav sample_sample.prm', ...
    'jrc3 export-spkwav sample_sample.prm 1', ...
    'jrc3 export-spkamp sample_sample.prm', ...
    'jrc3 export-spkamp sample_sample.prm 1', ...
    'jrc3 export-jrc1 sample_sample.prm', ...
    'jrc3 export-fet sample_sample.prm', ...
    'jrc3 plot-activity sample_sample.prm', ... %     'jrc3 kilosort sample_sample.prm', ...
    'jrc3 clear sample_sample.prm', ...
    'jrc3 traces-test sample_sample.prm', ...
    'jrc3 detectsort sample_sample.prm', ...
    'jrc3 auto sample_sample.prm', ...
    'jrc3 manual-test sample_sample.prm', ...
    }; %last one should be the manual test

    if ~isempty(vcArg1)
        switch lower(vcArg1)
            case {'show', 'info', 'list', 'help'}
            arrayfun(@(i)fprintf('%d: %s\n', i, csCmd{i}), 1:numel(csCmd));
            return;
            case {'manual', 'ui', 'ui-manual'}
            iTest = numel(csCmd); % + [-1,0];
            case {'traces', 'ui-traces'}
            iTest = numel(csCmd)-2; % second last
            otherwise
            iTest = str2num(vcArg1);
        end
        fprintf('Running test %s: %s\n', vcArg1, csCmd{iTest});
        csCmd = csCmd(iTest);
    end

    vlPass = false(size(csCmd));
    [csError, cS_prof] = deal(cell(size(csCmd)));
    vrRunTime = zeros(size(csCmd));
    for iCmd = 1:numel(csCmd)
        eval('close all; fprintf(''\n\n'');'); %clear memory
        fprintf('Test %d/%d: %s\n', iCmd, numel(csCmd), csCmd{iCmd});
        t1 = tic;
        profile('on');
        fDebug_ui = 1;
        setUserData(fDebug_ui);
        try
            if any(csCmd{iCmd} == '(' | csCmd{iCmd} == ';') %it's a function
                evalin('base', csCmd{iCmd}); %run profiler
            else % captured by profile
                csCmd1 = strsplit(csCmd{iCmd}, ' ');
                feval(csCmd1{:});
            end
            vlPass(iCmd) = 1; %passed test
        catch
            csError{iCmd} = lasterr();
            fprintf(2, '\tTest %d/%d failed\n', iCmd, numel(csCmd));
        end
        vrRunTime(iCmd) = toc(t1);
        cS_prof{iCmd} = profile('info');
    end
    nFailed = sum(~vlPass);

    fprintf('Unit test summary: %d/%d failed.\n', sum(~vlPass), numel(vlPass));
    for iCmd = 1:numel(csCmd)
        if vlPass(iCmd)
            fprintf('\tTest %d/%d (''%s'') took %0.1fs.\n', iCmd, numel(csCmd), csCmd{iCmd}, vrRunTime(iCmd));
        else
            fprintf(2, '\tTest %d/%d (''%s'') failed:%s\n', iCmd, numel(csCmd), csCmd{iCmd}, csError{iCmd});
        end
    end

    if numel(cS_prof)>1
        assignWorkspace_(cS_prof);
        disp('To view profile, run: profview(0, cS_prof{iTest});');
    else
        profview(0, cS_prof{1});
    end
    fDebug_ui = [];
    setUserData(fDebug_ui);
end %func
