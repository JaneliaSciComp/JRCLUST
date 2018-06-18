%--------------------------------------------------------------------------
function Sfile = importTSF_header_(arg1)
    % Sfile = importTSF_header(vcFname) %pass file name string
    % Sfile = importTSF_header(fid) %pass fid

    if ischar(arg1)
        fid = fopen(arg1, 'r');
        vcFname = arg1;
    else
        fid = arg1;
        vcFname = [];
    end

    % 16 bytes
    header = fread(fid, 16, 'char*1=>char');

    % 4x5 bytes
    iformat = fread(fid, 1, 'int32');
    SampleFrequency = fread(fid, 1, 'int32');
    n_electrodes = fread(fid, 1, 'int32');
    n_vd_samples = fread(fid, 1, 'int32');
    vscale_HP = fread(fid, 1, 'single');

    % electrode info
    Siteloc = zeros(2, n_electrodes, 'int16'); %2 x n_elec
    Readloc = zeros(1, n_electrodes, 'int32'); %read location
    for i_electrode = 1:n_electrodes
        Siteloc(:, i_electrode) = fread(fid, 2, 'int16');
        Readloc(i_electrode) = fread(fid, 1, 'int32');
    end
    if ~isempty(vcFname), fclose(fid); end

    Sfile = struct('header', header, 'iformat', iformat, ...
    'sRateHz', SampleFrequency, 'nChans', n_electrodes, ...
    'n_vd_samples', n_vd_samples, 'vscale_HP', vscale_HP, ...
    'Siteloc', Siteloc, 'tLoaded', n_vd_samples/SampleFrequency, 'Readloc', Readloc);
end %func
