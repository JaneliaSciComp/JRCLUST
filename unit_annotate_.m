%--------------------------------------------------------------------------
function unit_annotate_(hObject, event, vcLabel)
    S0 = get(0, 'UserData');
    S_clu = S0.S_clu;
    iClu1 = S0.primarySelectedCluster;
    if ~isfield(S_clu, 'clusterNotes'), S_clu.clusterNotes = cell(S_clu.nClusters, 1); end
    if nargin==3
        if isempty(vcLabel), vcLabel='';
        elseif vcLabel(1) == '='
            if ~isempty(S0.secondarySelectedCluster)
                vcLabel = sprintf('=%d', S0.secondarySelectedCluster);
            else
                msgbox_('Right-click another unit to set equal to.');
                return;
                vcLabel = '';
            end
        end
        S0.S_clu.clusterNotes{iClu1} = vcLabel;
    else
        vcNote1 = S_clu.clusterNotes{iClu1};
        if isempty(vcNote1), vcNote1=''; end
        csAns = inputdlg_(sprintf('Clu%d', iClu1), 'Annotation', 1, {vcNote1});
        if isempty(csAns), return; end
        vcLabel = csAns{1};
        S0.S_clu.clusterNotes{iClu1} = vcLabel;
    end

    % set(0, 'UserData', S0);
    clu_info_(S0); %update label
    save_log_(sprintf('annotate %d %s', iClu1, vcLabel), S0);
    % update cluster
end %func
