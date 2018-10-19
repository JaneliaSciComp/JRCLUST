%--------------------------------------------------------------------------
function export_csv_(varargin)
    % export_csv_(hObject, event)
    % if nargin<2,
    fZeroIndex = 0; %zero-based index export (disable to export as matlab index starting from 1)

    % S0 = get(0, 'UserData');
    if nargin==2
        [S0, P, S_clu] = get0_();
    elseif nargin==1
        P = varargin{1};
        vcFile_prm = P.vcFile_prm;
        S0 = load_cached_(P, 0);
        if isempty(S0), fprintf(2, 'Cannot find _jrc.mat.\n'); return; end %exit if file doesn't exist
        P = S0.P;
    end

    % vcFile_clu = subsFileExt(P.vcFile_prm, '_clu.mat');
    % Sclu = load(vcFile_clu); %load Sclu
    % if isfield(Sclu, 'Sclu'), Sclu = Sclu.Sclu; end

    if isfield(S0, 'S_clu')
        viClu = double(S0.S_clu.viClu);
    else
        fprintf(2, 'Cannot find S_clu.\n');
    end
    vrTime = double(S0.viTime_spk) / P.sRateHz;
    viSite = double(S0.viSite_spk) - fZeroIndex; %zero base

    vcFile_csv = subsFileExt_(P.vcFile_prm, '.csv');
    dlmwrite(vcFile_csv, [vrTime(:), viClu(:), viSite(:)], 'precision', 9);
    fprintf('Wrote to %s. Columns:\n', vcFile_csv);
    fprintf('\tColumn 1: Spike time (s)\n');
    fprintf('\tColumn 2: Unit# (positive #: valid units, 0: noise cluster, negative #: deleted clusters)\n');
    fprintf('\tColumn 3: Site# (starts with 1)\n');
end %func
