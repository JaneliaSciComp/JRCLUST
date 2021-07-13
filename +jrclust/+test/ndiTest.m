classdef ndiTest < matlab.uitest.TestCase & matlab.mock.TestCase
	methods(Test)
		function testNDIBootstrapAndDetect(tc)
			% use mock session
			S = ndi.session.mock();
			P = S.getprobes('type','n-trode');
			E = S.getelements('element.type','n-trode');
			E = E{1};
			jrc('bootstrap','ndi',S,E,'noShow');

			paramdir = dir([S.path() filesep '.JRCLUST' filesep '*_|_*'])

			paramfile = [S.path() filesep '.JRCLUST' filesep paramdir(1).name ...
				filesep 'jrclust.prm'];

			eval(['jrc detect ' paramfile]);
		end;
	end;
end
