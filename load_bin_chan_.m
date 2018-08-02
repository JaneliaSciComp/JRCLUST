%--------------------------------------------------------------------------
% 12/20/17 JJJ: support combined recordings and header containing files
% 12/15/17 JJJ: Load a single channel in memory efficient way
function vrWav = load_bin_chan_(P, iChan)
    % vrWav = load_bin_chan_(P, 0): return average of all valid channels
    % vrWav = load_bin_chan_(P, iChan): return specific channel
    if nargin<2, iChan = 0; end
    vrWav = [];

    % Join multiple recordings if P.multiFilenames is set
    if isempty(P.vcFile) && ~isempty(P.multiFilenames)
        csFile_bin = filter_files_(P.multiFilenames);
        cvrWav = cell(numel(csFile_bin), 1);
        P_ = setfield(P, 'multiFilenames', {});
        for iFile=1:numel(csFile_bin)
            P_ = setfield(P, 'vcFile', csFile_bin{iFile});
            cvrWav{iFile} = load_bin_chan_(P_, iChan);
            fprintf('.');
        end
        vrWav = cell2mat_(cvrWav);
        return;
    end

    try
        vrWav = load_bin_paged_(P, iChan);
    catch
        disperr_();
        return;
    end
end %func
