%--------------------------------------------------------------------------
function [vlIn_spk, mrFet, vhAx] = auto_split_wav_(mrSpkWav, mrFet, nSplits)
    % TODO: ask users number of clusters and split multi-way
    %Make automatic split of clusters using PCA + kmeans clustering
    %  input Sclu, trSpkWav and cluster_id of the cluster you want to cut

    if nargin<2, mrFet = []; end
    if nargin<3, nSplits = 2; end
    if isempty(mrFet)
        [~,pcaFet,vrD] = pca(double(mrSpkWav'), 'Centered', 1, 'NumComponents', 3);
        mrFet = pcaFet;
    else
        [~,pcaFet,vrD] = pca(double(mrSpkWav'), 'Centered', 1, 'NumComponents', 3);
        mrFet = double([mrFet pcaFet]);
    end

    inClust = MikeSplit(mrSpkWav, mrFet, nSplits);
    mrFet = pcaFet;
    nSpks = size(mrSpkWav,2);
    % nSplit = preview_split_(mrSpkWav1);
    % if isnan(nSplit), return; end

    spFig = resize_figure_([], [.5 0 .5 1]);
    figure(spFig);
    vhAx = zeros(4,1);
    for iAx=1:numel(vhAx)
        vhAx(iAx) = subplot(2,2,iAx); hold on;
    end
    % plot(vhAx(1), mrFet(:,1), mrFet(:,2), '.'); xylabel_(vhAx(1), 'PC1', 'PC2', 'PC1 vs PC2');
    % plot(vhAx(2), mrFet(:,3), mrFet(:,2), '.'); xylabel_(vhAx(2), 'PC3', 'PC2', 'PC3 vs PC2');
    % plot(vhAx(3), mrFet(:,1), mrFet(:,3), '.'); xylabel_(vhAx(3), 'PC1', 'PC3', 'PC1 vs PC3');
    drawnow;

    % Ask how many clusters there are
    try
        % kmean clustering into 2
        idx = kmeans(mrFet, nSplits);
        d12 = mad_dist_(mrFet(idx==1,:)', mrFet(idx==2,:)');
        fprintf('mad_dist: %f\n', d12);
        % idx = kmeans([pca_1,pca_2], NUM_SPLIT);
        % vlIn_spk = logical(idx-1);
        vlIn_spk = inClust;
    catch
        %         msgbox('Too few spikes to auto-split');
        vlIn_spk = false(nSpks,1);
        vlIn_spk(1:end/2) = true;
        d12 = nan;
        return;
    end
    vlOut_spk = ~vlIn_spk;
    plot(vhAx(1), mrFet(vlIn_spk,1), mrFet(vlIn_spk,2), 'b.', mrFet(vlOut_spk,1), mrFet(vlOut_spk,2), 'r.');
    plot(vhAx(2), mrFet(vlIn_spk,3), mrFet(vlIn_spk,2), 'b.', mrFet(vlOut_spk,3), mrFet(vlOut_spk,2), 'r.');
    plot(vhAx(3), mrFet(vlIn_spk,1), mrFet(vlIn_spk,3), 'b.', mrFet(vlOut_spk,1), mrFet(vlOut_spk,3), 'r.');

    min_y=min(reshape(mrSpkWav,1,[]));
    max_y=max(reshape(mrSpkWav,1,[]));
    viSpk1 = subsample_vr_(1:size(mrSpkWav,1), 1000); % show 1000 sites only

    hold(vhAx(4), 'on');
    title(vhAx(4), 'Mean spike waveforms');
    plot(vhAx(4), mean(mrSpkWav(viSpk1,vlIn_spk),2),'b');
    plot(vhAx(4), mean(mrSpkWav(viSpk1,vlOut_spk),2),'r');
    ylim_([min_y max_y]);
end %func
