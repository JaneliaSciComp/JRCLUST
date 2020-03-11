function [simScoreCorr,simScoreAmp,bestLag] = waveformSimScore(means,max_lag,sites)
% similarity scores using max cross-correlation (up to max_lag) between mean waveforms

% means is timepoints x channels x clusters
% sites can be used to index which set of channels is included for each
% cluster and is size channels x clusters. This is useful for high-density
% arrays where there is no point including all channels for each cluster.

[n_time_points,n_sites,n_clusters] = size(means);
if nargin<3
    sites = repmat([1:n_sites]',1,n_clusters);
else
    for i=1:n_clusters % sort them by sites, rather than maximum size
        [sites(:,i),sort_idx_i] = sort(sites(:,i));
        means(:,:,i) = means(:,sort_idx_i,i); 
    end
end
[simScoreCorr,simScoreAmp,bestLag] = deal(zeros(n_clusters,n_clusters));
L = (n_time_points*n_sites);
for i=1:n_clusters        
    if mod(i,10)==1
        tic;
    end   
    if nargin<3
        waveform_i = means(:,:,i);
        waveform_i = waveform_i(:);
        waveform_i = waveform_i - sum(waveform_i)./L;
        di = sqrt(sum(waveform_i.^2));
    end
    for j =1:n_clusters
        if i<j
            if nargin>2
                common_sites = intersect(sites(:,j),sites(:,i));
                if isempty(common_sites)
                    continue
                end
                waveform_j = means(:,ismember(sites(:,j),common_sites),j);
                waveform_i = means(:,ismember(sites(:,i),common_sites),i);   
                waveform_i = waveform_i(:);
                waveform_i = waveform_i - sum(waveform_i)./L;
                di = sqrt(sum(waveform_i.^2));                    
            else
               waveform_j = means(:,:,j); 
            end
            amp_i = sqrt(sum(waveform_i.^2));
            amp_j = sqrt(sum(waveform_j(:).^2));
            simScoreAmp(i,j) = (amp_i-amp_j)./(amp_i+amp_j);            
            for k = 1:(2*max_lag+1)
                    curr_j = zeros(size(waveform_j));
                    if max_lag==0
                        lag=0;
                    else
                        lag = k - max_lag -1;
                    end
                    if lag<0
                        curr_j(1+abs(lag):end,:) = waveform_j(1:end-abs(lag),:);
                    elseif lag>0
                        curr_j(1:end-abs(lag),:) = waveform_j(1+abs(lag):end,:);                    
                    else
                        curr_j = waveform_j;                    
                    end
                    curr_j = curr_j(:);
                    curr_j = curr_j - sum(curr_j)./L;    
                    coef = waveform_i' * curr_j;
                    dj = sqrt(sum(curr_j.^2));
                    coef = coef./(di*dj);
                    if coef>simScoreCorr(i,j)
                        simScoreCorr(i,j) = coef;
                        bestLag(i,j)=lag;
                    end
                    if coef==1
                        continue
                    end
             end
        end
    end
    if mod(i,10)==0
        fprintf('Finished clusters %g to %g in %s.\n',i-9,i,timestr(toc));
    end    
end
simScoreCorr = simScoreCorr+simScoreCorr' + eye(n_clusters); % symmetrize
simScoreAmp = simScoreAmp+simScoreAmp'; % symmetrize
bestLag = bestLag + bestLag';            


end