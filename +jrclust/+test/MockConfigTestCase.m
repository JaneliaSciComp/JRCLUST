classdef (Abstract) MockConfigTestCase < matlab.mock.TestCase
    %MOCKCONFIGTESTCASE Superclass of tests for all objects which require
    %a Config member.

    %% FIRST-CLASS PROPERTIES
    properties
        hCfg;                           % mock jrclust.Config object
        hCfgBehavior;                   % behavior object for mock Config
        histFile = tempname();          % temporary history file
        resFile = [tempname() '.mat'];  % temporary res file
        nSpikes = 8192;                 % 2^13 spikes
        nSites = 64;                    % 64 sites
    end
    
    %% SETUP METHODS
    methods (TestClassSetup)
        function setupConfig(obj)
            %SETUPCONFIG Create a mock jrclust.Config object with just the
            % necessary properties.
            import matlab.mock.actions.AssignOutputs;
            import matlab.mock.actions.Invoke;

            params = jrclust.utils.getDefaultParams();
            defaultParamNames = fieldnames(params)';

            paramNames = [defaultParamNames ...
                {'testRun', 'bytesPerSample', 'histFile', 'refracIntSamp', 'resFile', 'sessionName'}...
                {'evtWindowSamp', 'evtWindowRawSamp', 'nSites', 'nSitesEvt', 'siteNeighbors'}];

            [obj.hCfg, obj.hCfgBehavior] = obj.createMock( ...
                'AddedProperties', paramNames, ...
                'AddedMethods', ["getOr", ...
                                 "isa", ...
                                 "resetTemporaryParams", ...
                                 "setTemporaryParams", ...
                                 "updateLog", ...
                                ] ...
                );

            % set default param values
            for i = 1:numel(defaultParamNames)
                paramName = defaultParamNames{i};
                param = params.(paramName);

                val = param.default_value;
                if isfield(param.validation, 'postapply')
                    hFun = eval(param.validation.postapply);
                    val = hFun(val);
                end
                obj.assignOutputsWhen(get(obj.hCfgBehavior.(paramName)), val);
            end

            % bypass prompts
            obj.assignOutputsWhen(get(obj.hCfgBehavior.testRun), true);

            % set unavoidably user-defined values
            nSiteDir = 7;
            nSitesExcl = 2;
            
            obj.assignOutputsWhen(get(obj.hCfgBehavior.histFile), obj.histFile);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSiteDir), nSiteDir);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSitesExcl), nSitesExcl);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.resFile), obj.resFile);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.siteMap), (1:obj.nSites)');
            obj.assignOutputsWhen(get(obj.hCfgBehavior.siteLoc), rand(obj.nSites, 2));

            % set derived param values
            obj.assignOutputsWhen(get(obj.hCfgBehavior.bytesPerSample), 2); % default int16
            obj.assignOutputsWhen(get(obj.hCfgBehavior.evtWindowSamp), ...
                round(params.evtWindow.default_value * params.sampleRate.default_value / 1000));

            obj.assignOutputsWhen(get(obj.hCfgBehavior.evtWindowRawSamp), ...
                round(params.evtWindowRaw.default_value * params.sampleRate.default_value / 1000));

            obj.assignOutputsWhen(get(obj.hCfgBehavior.figPos), ...
                obj.defaultFigPos(params.figList.default_value));
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSites), obj.nSites);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.nSitesEvt), ...
                1 + 2*nSiteDir - nSitesExcl);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.refracIntSamp), ...
                round(params.refracInt.default_value * params.sampleRate.default_value / 1000));

            obj.assignOutputsWhen(get(obj.hCfgBehavior.sessionName), 'test');

            siteNeighbors = zeros(1 + 2*nSiteDir, obj.nSites);
            for i = 1:obj.nSites
                siteNeighbors(:, i) = mod((i-1:i + 2*nSiteDir - 1)', obj.nSites) + 1;
            end
            obj.assignOutputsWhen(get(obj.hCfgBehavior.siteNeighbors), siteNeighbors);

            % set method return values
            when(withAnyInputs(obj.hCfgBehavior.getOr), Invoke(@(varargin) obj.getOr(varargin)));
            when(withAnyInputs(obj.hCfgBehavior.isa), Invoke(@(varargin) obj.hCfgIsA(varargin)));
        end
        
        function setupProps(obj)
            % touch histFile
            fclose(fopen(obj.histFile, 'w'));

            % touch resFile
            fclose(fopen(obj.resFile, 'w'));
        end
    end
    
    %% TEARDOWN METHODS
    methods (TestClassTeardown)
        function rmHistFile(obj)
            fclose all;
            if exist(obj.histFile, 'file') == 2
                delete(obj.histFile);
            end
        end
    end
    
    %% MOCK METHODS
    methods (Access = protected)
        function val = getOr(varargin)
            val = [];
            args = varargin{2};
            if strcmp(args{2}, 'testRun')
                val = true;
            elseif numel(args) == 3
                val = args{3};
            end
        end
        
        function val = hCfgIsA(varargin)
            args = varargin{2};
            val = strcmp(args{end}, 'jrclust.Config');
        end
    end

    %% HELPER METHODS
    methods (Access = protected)
        function figPos = defaultFigPos(~, figList)
            figPos = cell(1, numel(figList));
            hasFigRD = ismember('FigRD', figList);

            for f = 1:length(figList)
               switch figList{f}
                   case 'FigCorr'
                       if hasFigRD
                           figPos{f} = [.85 .25 .15 .25];
                       else
                           figPos{f} = [.85 .2 .15 .27];
                       end

                   case 'FigHist'
                       if hasFigRD
                           figPos{f} = [.85 .75 .15 .25];
                       else
                           figPos{f} = [.85 .73 .15 .27];
                       end

                   case 'FigISI'
                       if hasFigRD
                           figPos{f} = [.85 .5 .15 .25];
                       else
                           figPos{f} = [.85 .47 .15 .26];
                       end

                   case 'FigMap'
                       figPos{f} = [0 .5 .15 .5];

                   case 'FigPos'
                       figPos{f} = [0 0 .15 .5];

                   case 'FigProj'
                       figPos{f} = [.5 .2 .35 .5];

                   case 'FigRD'
                       figPos{f} = [.85 0 .15 .25];

                   case 'FigSim'
                       figPos{f} = [.5 .7 .35 .3];

                   case 'FigTime'
                       if hasFigRD
                           figPos{f} = [.15 0 .7 .2];
                       else
                           figPos{f} = [.15 0 .85 .2];
                       end

                   case 'FigWav'
                       figPos{f} = [.15 .2 .35 .8];
               end
            end
        end
    end
end

