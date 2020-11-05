function assign (varargin)
% assign - make a list of assignments (matlab 5 or higher)
%
%	ASSIGN('VAR1', VAL1, 'VAR2', VAL2, ...) makes the assignments 
%	VAR1 = VAL1; VAR2 = VAL2; ... in the caller's workspace.
%
%	This is most useful when passing in an option list to a
%	function.  Thus in the function which starts:
%		function foo(x,y,varargin)
%		z = 0;
%		assign(varargin{:});
%	the variable z can be given a non-default value by calling the
%	function like so: foo(x,y,'z',4);
%
%       If the input is a single structure, then the structure is converted
%       to a set of NAME/VALUE pairs and interpreted as 'VAR1', VAL1, etc,
%       where VAR1 is the first field name of the input, VAL1 is the value of the field name,
%       etc.


if numel(varargin)==1,
    if isstruct(varargin{1}),
        varargin = struct2namevaluepair(varargin{1});
    elseif iscell(varargin{1}),
        varargin = varargin{1};
    end;
end;

vars = {varargin{1:2:end}};
vals = {varargin{2:2:end}};

% use length(vals) not length(vars) in case of odd number of arguments
for i = 1:length(vals), % changed 1 to a 2
  assignin('caller', vars{i}, vals{i});
end
