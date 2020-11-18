classdef SaveStructTest < matlab.unittest.TestCase
    %SAVESTRUCTTEST Tests for jrclust.utils.saveStruct.
    properties
        testStruct;
    end

    %% SETUP METHODS
    methods (TestClassSetup)
        function setupStructData(obj)
            obj.testStruct = struct('foo', 'bar', 'baz', 1, 'qux', {'quux'});
        end
    end

    %% TEST METHODS
    methods (Test)
        function saveSimpleOK(obj)
            %SAVESIMPLEOK Test that saving a simple struct to a simple path
            %works.
            filename = tempname();
            obj.assertTrue(jrclust.utils.saveStruct(obj.testStruct, filename));

            S = load(filename);
            obj.assertEqual(S, obj.testStruct);
        end

        function catchLongFilename(obj)
            %CATCHLONGFILENAME Test that trying to save a struct to a very
            %long path (on Windows) will save to a temp file and copy it
            %over.
            if ispc
                % Windows' limit without registry hacking is an absurd 260
                % characters
                filelen = 260 - numel(tempdir());
                longfilename = fullfile(tempdir(), [repmat('a', 1, filelen) '.mat']);
                obj.assertTrue(jrclust.utils.saveStruct(obj.testStruct, longfilename));

                S = load(['\\?\' longfilename], '-mat');
                obj.assertEqual(S, obj.testStruct);
            else
                obj.assertTrue(true);
            end
        end
    end
end

