%--------------------------------------------------------------------------
function [mrWav, S_tsf] = importTSF_(fname, varargin)
    fid = fopen(fname, 'r');
    S_tsf = importTSF_header_(fid);
    n_vd_samples = S_tsf.n_vd_samples;
    n_electrodes = S_tsf.nChans;
    mrWav = reshape(fread(fid, n_vd_samples * n_electrodes, '*int16'), [n_vd_samples,n_electrodes]);
    fclose(fid);
end % function
