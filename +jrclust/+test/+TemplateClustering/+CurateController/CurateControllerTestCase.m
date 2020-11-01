classdef (Abstract) CurateControllerTestCase < jrclust.test.TemplateClustering.TemplateClusteringTestCase
    %CURATECONTROLLERTESTCASE Superclass of tests for CurateController
    %objects, or objects which require a CurateController member.
    %   Uses TemplateClustering for hClust object.

    %% FIRST-CLASS PROPS
    properties
        hCurate;    % CurateController object
    end

    %% DEPENDENT PROPS
    properties (Dependent)
        selected;   % CurateController's selected unit(s).
    end
    
    %% SETUP METHODS
    methods (TestClassSetup)
        function setupConfig(obj)
            %SETUPCONFIG CurateController does not use FigRD for
            %TemplateClustering, so remove it from figList.
            setupConfig@jrclust.test.MockConfigTestCase(obj);

            import matlab.mock.actions.AssignOutputs;

            % skip FigRD in Curate tests
            figList = setdiff(obj.hCfg.figList, {'FigRD'});
            figPos = obj.defaultFigPos(figList);

            obj.assignOutputsWhen(get(obj.hCfgBehavior.figList), figList);
            obj.assignOutputsWhen(get(obj.hCfgBehavior.figPos), figPos);
        end

        function setupProps(obj)
            %SETUPPROPS Create the necessary data for a mock Clustering and
            %CurateController.
            setupProps@jrclust.test.TemplateClustering.TemplateClusteringTestCase(obj);

            res = jrclust.utils.mergeStructs(obj.dRes, obj.sRes);
            res.hClust = obj.hClust;
            obj.hCurate = jrclust.curate.CurateController(res);
            obj.hCurate.beginSession();
        end
    end

    %% TEARDOWN METHODS
    methods (TestClassTeardown)
        function endSession(obj)
            obj.hCurate.endSession();
        end
    end

    methods (TestMethodTeardown)
        function reset(obj)
            %RESET Restore the Clustering and CurateController to their
            %initial states.
            reset@jrclust.test.TemplateClustering.TemplateClusteringTestCase(obj);

            % reset is called at test class setup, before hCurate is
            % instantiated
            if ~isempty(obj.hCurate)
                obj.hCurate.updateSelect(1);
                obj.hCurate.updateFigWav();
                obj.hCurate.updateFigRD();
                obj.hCurate.updateFigSim();
                obj.hCurate.updateHistMenu();
            end
        end
    end

    %% HELPER METHODS
    methods
        function updateSelect(obj, unitIds)
            %UPDATESELECT Update the CurateController's selected units.
            obj.hCurate.updateSelect(unitIds);
        end
    end

    %% GETTERS/SETTERS
    methods
        function s = get.selected(obj)
            s = obj.hCurate.selected;
        end
    end
end
