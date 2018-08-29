%--------------------------------------------------------------------------
function delete_annotate() %TW
% SNR based delete functionality
% Ask SNR
S0 = get(0, 'UserData');
[S_clu, P] = get0_('S_clu', 'P');

% Auto delete
figure_wait_(1); drawnow;

S_clu = delete_clu_(S_clu, find(strcmp(S_clu.clusterNotes, 'to_delete')));
setUserData(S_clu);
S0 = gui_update_();
figure_wait_(0);


save_log_(sprintf('delete-anotate'), S0);
end % function
