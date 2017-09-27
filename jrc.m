function varargout = jrc(varargin)
% calls jrclust

%jrclust(varargin{:});
fprintf('Running ''%s%sjrc3.m''\n', fileparts(mfilename('fullpath')), filesep());
jrc3(varargin{:});