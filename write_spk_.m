%--------------------------------------------------------------------------
% 10/11/17 JJJ: created and tested
function write_spk_(varargin) % DEPRECATED
    % [Usage]
    % write_spk_() % close and clear
    % write_spk_(vcFile_prm) %open file
    % write_spk_(spikeTraces, spikeWaveforms, spikeFeatures)
    persistent fid_raw fid_spk fid_fet

    switch nargin
        case 0
                fid_raw = fclose_(fid_raw);
                fid_spk = fclose_(fid_spk);
                fid_fet = fclose_(fid_fet);
        case 1
                vcFile_prm = varargin{1};
                fid_raw = fopen(strrep(vcFile_prm, '.prm', '_traces.bin'), 'W');
                fid_spk = fopen(strrep(vcFile_prm, '.prm', '_waveforms.bin'), 'W');
                fid_fet = fopen(strrep(vcFile_prm, '.prm', '_features.bin'), 'W');
        case 3
                fwrite_(fid_raw, varargin{1});
                fwrite_(fid_spk, varargin{2});
                fwrite_(fid_fet, varargin{3});
        otherwise
                disperr_('write_spk_: invalid nargin');
    end %switch
end % function
