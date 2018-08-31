%--------------------------------------------------------------------------
% 9/27/17 JJJ: Imported from SPARC grant
function [vrWav_aux, vrTime_aux] = load_aux_(P)

    [vrWav_aux, vrTime_aux] = deal([]);
    if ~isfield(P, 'vcFile') || isempty(P.vcFile), msgbox_('Multi-file mode is currently not supported'); return; end
    [~,~,vcExt] = fileparts(P.vcFile);
    vcFile_aux = getOr(P, 'vcFile_aux', '');
    if isempty(vcFile_aux)
        switch lower(vcExt)
            case '.ns5', vcFile_aux = subsFileExt_(P.vcFile, '.ns2');
            case {'.bin', '.dat'}, vcFile_aux = P.vcFile;
            otherwise
            fprintf(2, 'Unable to determine the aux file. You must manually specify "vcFile_aux".\n');
            return;
        end
    end
    if ~fileExists(vcFile_aux), return; end

    [~,~,vcExt_aux] = fileparts(vcFile_aux);
    switch lower(vcExt_aux)
        case '.ns2'
        iChan_aux = getOr(P, 'iChan_aux', 1);
        [mnWav_aux, hFile_aux, S_aux] = load_nsx_(vcFile_aux);
        scale_aux = hFile_aux.Entity(iChan_aux).Scale * P.vrScale_aux;
        vrWav_aux = single(mnWav_aux(iChan_aux,:)') * scale_aux;
        sampleRateHz_aux = S_aux.sampleRateHz;
        case '.mat'
        S_aux = load(vcFile_aux);
        csField_aux = fieldnames(S_aux);
        vrWav_aux = S_aux.(csField_aux{1});
        sampleRateHz_aux = getOr(P, 'sampleRateHz_aux', P.sampleRateHz_rate);
        case {'.dat', '.bin'}
        iChan_aux = getOr(P, 'iChan_aux', []);
        if isempty(iChan_aux), return; end
        mnWav_aux = load_bin_(vcFile_aux, P.dataType);
        vrWav_aux = single(mnWav_aux(iChan_aux:P.nChans:end)') * P.uV_per_bit * P.vrScale_aux;
        sampleRateHz_aux = getOr(P, 'sampleRateHz_aux', P.sampleRateHz);
        otherwise
        fprintf(2, 'vcFile_aux: unsupported file format: %s\n', vcExt_aux);
        return;
    end %switch
    if nargout>=2, vrTime_aux = single(1:numel(vrWav_aux))' / sampleRateHz_aux; end
end % function
