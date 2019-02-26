classdef SortTest < matlab.unittest.TestCase
    %SORTTEST Test spike clustering
    properties
        hJRC_CPU;
        hJRC_GPU;
    end

    properties (Dependent)
        res_CPU;
        res_GPU;
        hCfg_CPU;
        hCfg_GPU;
    end

    methods (TestClassSetup)
        function setupProps(obj)
            hCfg_CPU_ = jrclust.Config(fullfile(getenv('JRCTESTDATA'), 'large', 'test_large.prm'));
            obj.hJRC_CPU = JRC(hCfg_CPU_);
            obj.hCfg_CPU.testRun = 1;

            hCfg_GPU_ = jrclust.Config(fullfile(getenv('JRCTESTDATA'), 'large', 'test_large.prm'));
            obj.hJRC_GPU = JRC(hCfg_GPU_);
            obj.hCfg_GPU.testRun = 1;
        end

        function sortCPU(obj)
            obj.hCfg_CPU.useGPU = 0;
            obj.hJRC_CPU.sort();
        end

        function sortGPU(obj)
            obj.hCfg_GPU.useGPU = 1;
            obj.hJRC_GPU.sort();
        end
    end

    methods (Test)
        function valsAgree(obj)
            %VALSAGREE Check CPU/GPU rho-delta computation is the same
            % rho should be the same
            obj.assertEqual(obj.res_CPU.spikeRho, obj.res_GPU.spikeRho, 'AbsTol', 1e-5);

            % reasonable people can disagree on very large values of delta
            deltaCPU = obj.res_CPU.spikeDelta;
            deltaGPU = obj.res_GPU.spikeDelta;
            obj.assertEqual(deltaCPU(deltaCPU < 1e16), deltaGPU(deltaGPU < 1e16), 'AbsTol', 1e-5);

            % nearest neighbors should be the same
            obj.assertEqual(obj.res_CPU.spikeNeigh, obj.res_GPU.spikeNeigh);

            % in the end the clustering should come out okay anyway
            obj.assertEqual(obj.res_CPU.spikeClusters, obj.res_GPU.spikeClusters);
        end

        function nonnegNoGaps(obj)
            %NONNEGNOGAPS Check clusters nonnegative and no gaps
            obj.assertGreaterThanOrEqual(obj.res_CPU.spikeClusters, 0);
            clusters = unique(obj.res_CPU.spikeClusters);

            % should be 1:nClusters with no gaps
            obj.assertEqual(numel(clusters(clusters > 0)), max(clusters));
        end
    end

    %% GETTERS/SETTERS
    methods
        function hCfg_CPU = get.hCfg_CPU(obj)
            hCfg_CPU = obj.hJRC_CPU.hCfg;
        end
        function set.hCfg_CPU(obj, hCfg)
            obj.hJRC_CPU.hCfg = hCfg;
        end

        function hCfg_GPU = get.hCfg_GPU(obj)
            hCfg_GPU = obj.hJRC_GPU.hCfg;
        end
        function set.hCfg_GPU(obj, hCfg)
            obj.hJRC_GPU.hCfg = hCfg;
        end

        function res_CPU = get.res_CPU(obj)
            res_CPU = obj.hJRC_CPU.res;
        end

        function res_GPU = get.res_GPU(obj)
            res_GPU = obj.hJRC_GPU.res;
        end
    end
end


