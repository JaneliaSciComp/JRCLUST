classdef DeprecationTest < matlab.unittest.TestCase
    %DEPRECATIONTESTS test that deprecated commands are indeed deprecated
    
    methods (Test)
        function compileKsort(testCase)
            hJRC = jrc('compile-ksort');
            testCase.assertFalse(hJRC.doRun);
        end

        function doc(testCase)
            hJRC = jrc('doc');
            testCase.assertFalse(hJRC.doRun);
        end

        function docEdit(testCase)
            hJRC = jrc('doc-edit');
            testCase.assertFalse(hJRC.doRun);
        end

        function download(testCase)
            hJRC = jrc('download');
            testCase.assertFalse(hJRC.doRun);
        end

        function gitPull(testCase)
            hJRC = jrc('git-pull');
            testCase.assertFalse(hJRC.doRun);
        end

        function importKilosortSort(testCase)
            hJRC = jrc('import-kilosort-sort');
            testCase.assertFalse(hJRC.doRun);
        end

        function importKsortSort(testCase)
            hJRC = jrc('import-ksort-sort');
            testCase.assertFalse(hJRC.doRun);
        end

        function install(testCase)
            hJRC = jrc('install');
            testCase.assertFalse(hJRC.doRun);
        end

        function issue(testCase)
            hJRC = jrc('issue');
            testCase.assertFalse(hJRC.doRun);
        end

        function kilosort(testCase)
            hJRC = jrc('kilosort');
            testCase.assertFalse(hJRC.doRun);
        end

        function kilosortVerify(testCase)
            hJRC = jrc('kilosort-verify');
            testCase.assertFalse(hJRC.doRun);
        end
        
        function ksort(testCase)
            hJRC = jrc('ksort');
            testCase.assertFalse(hJRC.doRun);
        end

        function ksortVerify(testCase)
            hJRC = jrc('ksort-verify');
            testCase.assertFalse(hJRC.doRun);
        end

        function prmSet(testCase)
            hJRC = jrc('set');
            testCase.assertFalse(hJRC.doRun);
        end

        function setPrm(testCase)
            hJRC = jrc('set-prm');
            testCase.assertFalse(hJRC.doRun);
        end

        function setPrm2(testCase)
            hJRC = jrc('setprm');
            testCase.assertFalse(hJRC.doRun);
        end

        function update(testCase)
            hJRC = jrc('update');
            testCase.assertFalse(hJRC.doRun);
        end

        function which(testCase)
            hJRC = jrc('which');
            testCase.assertFalse(hJRC.doRun);
        end

        function wiki(testCase)
            hJRC = jrc('wiki');
            testCase.assertFalse(hJRC.doRun);
        end

        function wikiDownload(testCase)
            hJRC = jrc('wiki-download');
            testCase.assertFalse(hJRC.doRun);
        end
    end
end

