classdef (Abstract) CurateControllerTestCase < jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase
    %CURATECONTROLLERTESTCASE Superclass of tests for CurateController
    %objects, or objects which require a CurateController member.
    %   Uses DensityPeakClustering for hClust object.

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
        function setupProps(obj)
            %SETUPPROPS Create the necessary data for a mock Clustering and
            %CurateController.
            setupProps@jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase(obj);

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
            reset@jrclust.test.DensityPeakClustering.DensityPeakClusteringTestCase(obj);

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
