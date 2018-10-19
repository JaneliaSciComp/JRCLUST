%--------------------------------------------------------------------------
function nFailed = unit_test_(vcArg)
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

    if nargin < 1
        vcArg = '';
    end

    basedir = fileparts(fileparts(mfilename('fullpath')));
    cd(fullfile(basedir, 'test'));

    if ~exist_file_('sample.bin') || ~exist_file_('sample.meta')
        fprintf(2, 'test data missing');
        fprintf(2, 'find test data at https://drive.google.com/drive/folders/1-UTasZWB0TwFFFV49jSrpRPHmtve34O0?usp=sharing');
        nFailed = 0;
        return;
    end

    nFailed = 0;
    profile('clear'); %reset profile stats
    csCmd = {...
    'close all; clear all;', ... %start from blank
    'jrc clear sample_sample.prm', ...
    'jrc compile', ...
    'jrc probe sample.prb', ...
    'jrc makeprm sample_list.txt sample.prb', ...
    'jrc makeprm sample.bin sample.prb', ...
    'jrc probe sample_sample.prm', ...
    'jrc import-lfp sample_sample.prm', ...
    'jrc traces-lfp sample_sample.prm', ...
    'jrc traces sample_sample.prm', 'jrc traces', ...
    'jrc traces-test sample_sample.prm', ...
    'jrc preview-test sample_sample.prm', ...
    'jrc makeprm sample.bin sample.prb', ...
    'jrc detectsort sample_sample.prm', 'jrc clear sample_sample.prm', ...
    'jrc probe sample_sample_merge.prm', 'jrc detectsort sample_sample_merge.prm', 'jrc clear sample_sample_merge.prm', ... %multishank, multifile test
    'jrc detect sample_sample.prm', 'jrc sort sample_sample.prm', ...
    'jrc export-csv sample_sample.prm', ...
    'jrc export-quality sample_sample.prm', ...
    'jrc export-spkwav sample_sample.prm', ...
    'jrc export-spkwav sample_sample.prm 1', ...
    'jrc export-spkamp sample_sample.prm', ...
    'jrc export-spkamp sample_sample.prm 1', ...
    'jrc export-jrc1 sample_sample.prm', ...
    'jrc export-fet sample_sample.prm', ...
    'jrc plot-activity sample_sample.prm', ... %     'jrc kilosort sample_sample.prm', ...
    'jrc clear sample_sample.prm', ...
    'jrc traces-test sample_sample.prm', ...
    'jrc detectsort sample_sample.prm', ...
    'jrc auto sample_sample.prm', ...
    'jrc manual-test sample_sample.prm', ...
    }; %last one should be the manual test

    if ~isempty(vcArg)
        switch lower(vcArg)
            case {'show', 'info', 'list', 'help'}
                arrayfun(@(i)fprintf('%d: %s\n', i, csCmd{i}), 1:numel(csCmd));
                return;

            case {'manual', 'ui', 'ui-manual'}
                iTest = numel(csCmd); % + [-1,0];

            case {'traces', 'ui-traces'}
                iTest = numel(csCmd)-2; % second last

            otherwise
                iTest = str2num(vcArg);
        end
        fprintf('Running test %s: %s\n', vcArg, csCmd{iTest});
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
        set0_(fDebug_ui);
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
    set0_(fDebug_ui);
end %func
