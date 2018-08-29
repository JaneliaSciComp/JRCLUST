%--------------------------------------------------------------------------
function S = imec3_imroTbl_(cSmeta)
    % Smeta has imroTbl

    vcDir_probe = 'C:\Dropbox (HHMI)\IMEC\SpikeGLX_Probe_Cal_Data\'; %this may crash. probe calibaration folder

    if isstruct(cSmeta), cSmeta = {cSmeta}; end %turn it into cell of struct
    % parse imroTbl
    cs_imroTbl = cellfun(@(S)S.imroTbl, cSmeta, 'UniformOutput', 0);
    cvn_imroTbl = cellfun(@(vc)textscan(vc, '%d', 'Delimiter', '( ),'), cs_imroTbl, 'UniformOutput', 0);
    cvn_imroTbl = cellfun(@(c)c{1}, cvn_imroTbl, 'UniformOutput', 0);
    S.viBank = cellfun(@(vn)vn(7), cvn_imroTbl);
    S.viRef = cellfun(@(vn)vn(8), cvn_imroTbl);
    S.vrGain_ap = single(cellfun(@(vn)vn(9), cvn_imroTbl));
    S.vrGain_lf = single(cellfun(@(vn)vn(10), cvn_imroTbl));
    S.nSites_bank = cvn_imroTbl{1}(4);

    Smeta1 = cSmeta{1};

    % correct gain
    nFiles = numel(S.viBank);
    nSites = numel(Smeta1.viSites);
    [mrScale_ap, mrScale_lf] = deal(ones(nSites, nFiles));
    S.vcProbeSN = sprintf('1%d%d', Smeta1.imProbeSN, Smeta1.imProbeOpt);
    % read gain correction
    vrGainCorr = ones(1, S.nSites_bank*4);
    if Smeta1.imProbeOpt ~= 2
        try
            vcFile_csv = sprintf('%s1%d%d\\Gain correction.csv', vcDir_probe, Smeta1.imProbeSN, Smeta1.imProbeOpt);
            try
                vrGainCorr = csvread(vcFile_csv, 1, 1);
            catch
                vrGainCorr = csvread(vcFile_csv, 0, 0);
            end
        catch
            ;
        end
    end

    % build scale
    for iFile = 1:nFiles
        vrGainCorr1 = vrGainCorr(Smeta1.viSites + double(S.nSites_bank*S.viBank(iFile)));
        mrScale_ap(:,iFile) = 1.2 * 1e6 / 2^10 / S.vrGain_ap(iFile) .* vrGainCorr1;
        mrScale_lf(:,iFile) = 1.2 * 1e6 / 2^10 / S.vrGain_lf(iFile) .* vrGainCorr1;
    end
    S.mrScale_ap = mrScale_ap;
    S.mrScale_lf = mrScale_lf;
end % function
