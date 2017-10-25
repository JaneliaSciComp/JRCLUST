function jrc(varargin)
% calls jrclust
warning off;
%jrclust(varargin{:});
fprintf('Running ''%s%sjrc3.m''\n', fileparts(mfilename('fullpath')), filesep());
jrc3(varargin{:});