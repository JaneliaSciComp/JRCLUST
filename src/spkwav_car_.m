%--------------------------------------------------------------------------
function [trWav2, mrWav_ref] = spkwav_car_(trWav2, P, nSites_spk, viSite2_spk)
    %function [trWav2, mrWav_ref] = trWav_car_sort_(trWav2, P)
    %  trWav2: nT x nSpk x nSites, single
    if nargin<3, nSites_spk = []; end
    if nargin<4, viSite2_spk = []; end
    vcSpkRef = get_set_(P, 'vcSpkRef', 'nmean');
    if strcmpi(vcSpkRef, 'nmean')
        fSort_car = 1;
    else
        fSort_car = -1; % no local subtraction
    end
    % fSort_car = get_set_(P, 'fSort_car', 1);
    if isempty(nSites_spk)
        nSites_spk = 1 + 2 * P.maxSite - P.nSites_ref; % size(tnWav_spk, 2);
    end
    if nSites_spk==1, mrWav_ref=[]; return; end

    % if P.nSites_ref==0, fSort_car = -1; end
    % nSites_spk = 1 + 2 * P.maxSite; % size(tnWav_spk, 2);
    dimm1 = size(trWav2);
    [nT_spk, nSpk] = deal(dimm1(1), dimm1(2));
    switch fSort_car
        case 1 % use n sites having the least SD as reference sites
        if isempty(viSite2_spk)
            %             viSite_ref_ = 2:nSites_spk;
            %             viSite_ref_ = ceil(nSites_spk/2):nSites_spk;
            viSite_ref_ = ceil(size(trWav2,3)/2):size(trWav2,3);
            mrWav_ref = mean(trWav2(:,:,viSite_ref_), 3);
        else
            trWav3 = trWav2(:,:,1:nSites_spk);
            trWav3(:,:,1) = 0;
            for iSpk1 = 1:numel(viSite2_spk)
                trWav3(:,iSpk1,viSite2_spk(iSpk1)) = 0;
            end
            mrWav_ref = sum(trWav3, 3) / (nSites_spk-2);
        end
        case 3
        mrWav_ref = mean(trWav2(:,:,2:end), 3);
        case 2
        [~, miSites_ref] = sort(squeeze_(var(trWav2)), 2, 'descend'); % use lest activities for ref
        %miSites_ref = miSites_ref(:,nSites_spk+1:end);
        miSites_ref = miSites_ref(:,3:nSites_spk);
        ti_dimm1 = repmat((1:nT_spk)', [1, nSpk, P.nSites_ref]);
        ti_dimm2 = repmat(1:nSpk, [nT_spk,1,P.nSites_ref]);
        ti_dimm3 = repmat(shiftdim(miSites_ref,-1), [nT_spk,1,1]);
        mrWav_ref = mean(trWav2(sub2ind(size(trWav2), ti_dimm1, ti_dimm2, ti_dimm3)),3);
        case 0
        mrWav_ref = mean(trWav2(:,:,nSites_spk+1:end), 3);
        case 4
        mrWav_ref = mean(trWav2(:,:,1:nSites_spk), 3);
        case 5
        mrWav_ref = mean(trWav2(:,:,2:end), 3);
        case -1
        mrWav_ref = [];
    end
    trWav2 = trWav2(:,:,1:nSites_spk);
    dimm2 = size(trWav2);
    if ismatrix(trWav2), dimm2(end+1) = 1; end
    if ~isempty(mrWav_ref)
        trWav2 = meanSubt_(reshape(bsxfun(@minus, reshape(trWav2,[],dimm2(3)), mrWav_ref(:)), dimm2));
        %     trWav2 = reshape(bsxfun(@minus, reshape(trWav2,[],dimm2(3)), mrWav_ref(:)), dimm2); % no meanSubt_
    else
        trWav2 = meanSubt_(trWav2);
    end
end %func
