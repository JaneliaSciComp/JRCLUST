% --------------------------------------------------------------------------
% function [S_clu, nMerges_clu] = S_clu_isosplit_(S_clu)
% % Jeremy Maglund inspired
% error('not implemented yet');
% % simple merge only, no cluster refinement
% % average cluster waveform based merging only, no individual spike split @TODO
% 
% % global tnWav_spk
% 
% nRepeat = 3;
% 
% P = get0_('P');
% % Build a distance matrix and identify the clostest pairs
% % comparisons_made = ;
% % final_pass = false;
% % labels = S_clu.viClu;
% % Kmax = S_clu.nClu;
% mrFet_clu = compute_centers_(S_clu); %nDim x nClu
% mrDist_clu = make_dists_matrix_(mrFet_clu, P); % Kmax x Kmax, inf diagonal
% 
% % preserve cluster index while merging and compact it later
% for iRepeat = 1:nRepeat
% %     active_labels = unique(labels);
% %     mrFet_mean_clu = fet_mean_clu_(S_clu);
% %     mrDist_clu = dist_clu_(S_clu, mrFet_mean_clu); % 
%     [inds1, inds2] = get_pairs_to_compare_(mrDist_clu);
%     if isempty(inds1), break; end
% 
%     [S_clu, mrFet_clu, vi_inds_merged] = compare_pairs_(S_clu, mrFet_clu, inds1, inds2);    
%     if isempty(vi_inds_merged), return; end %@TODO: may have to repartition later
%     
%     % update changed clusters and distance compared    
%     mrDist_clu = update_dists_matrix_(mrDist_clu, mrFet_clu, inds1, inds2, vi_inds_merged);    
% end %for
% 
% % [miClu_pairs, mrDist_clu] = clu_dist();
% % S_iso = struct('isocut_threshold', 1, 'refine_clusters', 1, 'max_iterations_per_pass', 500, 'whiten_cluster_pairs', 1, 'prevent_merge', false);
% % 
% % mrFet_center_clu = compute_centers_(mrFet, S_clu.viClu); 
% % final_pass = false;
% % comparisons_made=zeros(Kmax,Kmax);
% % while 1
% %     something_merged = false;
% %     clusters_changed_vec_in_pass=zeros(1,Kmax);
% %     for iLoop = 1:S_iso.max_iterations_per_pass
% %         
% %     end %for
% %     % zero out the comparisons made matrix only for those that have changed
% %     clusters_changed=find(clusters_changed_vec_in_pass);
% %     for j=1:length(clusters_changed)
% %         comparisons_made(clusters_changed(j),:)=0;
% %         comparisons_made(:,clusters_changed(j))=0;
% %     end;
% %     
% %     if (something_merged), final_pass=false; end;
% %     if (final_pass), break; end; % This was the final pass and nothing has merged
% %     if (~something_merged), final_pass=true; end; % If we are done, do one last pass for final redistributes
% % end
% end %func


%--------------------------------------------------------------------------
% function centers = compute_centers_(X,labels)
% % JJJ optimized, isosplit
% [M,N]=size(X);
% centers=zeros(M,N);
% counts=accumarray(labels',1,[N,1])';
% for m=1:M
%     centers(m,:)=accumarray(labels',X(m,:)',[N,1])';
% end
% viLabel = find(counts);
% centers(:,viLabel)= bsxfun(@rdivide, centers(:,viLabel), counts(viLabel));
% end %func


%--------------------------------------------------------------------------
% function mrFet_clu = compute_centers_(S_clu)
% % calculate mean cluster waveform, using unique waveforms
% global tnWav_spk
% dimm_spk = size(tnWav_spk);
% mrFet_clu = zeros(dimm_spk(2), S_clu.nClu, 'single');
% % my feature is average cluster waveform
% for iClu=1:S_clu.nClu
%     mrFet_clu(:,iClu) = ;
% end %for
% end %func


%--------------------------------------------------------------------------
% function mrDist_clu = make_dists_matrix_(mrFet_clu, P)
% 
% end %func


%--------------------------------------------------------------------------
% function [inds1, inds2] = get_pairs_to_compare_(mrDist_clu)
% % find mutual min-dist pairs, remove inf
% 
% end %func


%--------------------------------------------------------------------------
% function [S_clu, centers, viPair_merged] = compare_pairs_(S_clu, centers, viClu1, viClu2)
% % find mutual min-dist pairs: correlation based merging only
% nPairs = numel(viClu1);
% for iPair = 1:nPairs
%     
% end 
% 
% centers = update_centers_(centers);
% end %func


%--------------------------------------------------------------------------
% function mrDist_clu = update_dists_matrix_(mrDist_clu, inds1, inds2, vi_inds_merged)
% % find mutual min-dist pairs
% mrDist_clu(sub2ind(size(mrDist_clu), inds1, inds2)) = inf; % compared dist ignored
% 
% end %func


%--------------------------------------------------------------------------
% function [inds1,inds2] = get_pairs_to_compare_(dists, thresh)
% [~,N]=size(dists); % square matrix, nClu
% inds1=[];
% inds2=[];
% % dists=make_dists_matrix(centers);
% % dists(find(comparisons_made(:)))=inf;
% % for j=1:N
% %     dists(j,j)=inf;
% % end;
% % important to only take the mutal closest pairs -- unlike how we originally did it
% %something_changed=1;
% %while (something_changed)
%     %something_changed=0;
%     [~,best_inds]=min(dists,[],1);
%     for j=1:N
%         if (best_inds(j)>j)
%             if (best_inds(best_inds(j))==j) % mutual
%                 if (dists(j,best_inds(j))<thresh)
%                     inds1(end+1)=j;
%                     inds2(end+1)=best_inds(j);
%                     dists(j,:)=inf;
%                     dists(:,j)=inf;
%                     dists(best_inds(j),:)=inf;
%                     dists(:,best_inds(j))=inf;
%                     %something_changed=1;
%                 end
%             end
%         end     
%     end
% %end;
% end %func

