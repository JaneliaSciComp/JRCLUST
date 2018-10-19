%--------------------------------------------------------------------------
function batch_verify_(vcFile_batch, vcCommand)
    % batch process parameter files (.batch) file
    % Example
    %   jrc batch-verify skip my.batch
    %       just does the verification plot for all files in .batch file
    if ~exist(vcFile_batch, 'file'), fprintf(2, 'File does not exist\n'); return; end
    edit(vcFile_batch); %show script
    if nargin<2, vcCommand=[]; end
    if isempty(vcCommand), vcCommand='spikesort'; end
    csFiles_prm = load_batch_(vcFile_batch);

    if ~strcmpi(vcCommand, 'skip')
        for iFile=1:numel(csFiles_prm)
            try
                vcFile_prm1 = csFiles_prm{iFile};
                jrc('clear');
                jrc(vcCommand, vcFile_prm1);
                if isempty(strfind(vcCommand, 'verify'))
                    jrc('verify', vcFile_prm1); % try silent verify and collect result
                end
            catch
                disperr_();
            end
        end %for
    end
    fprintf('\nSummary for %s\n', vcFile_batch);
    % Collect data
    [cvrSnr, cvrFp, cvrFn, cvrAccuracy, cvnSite, cvnSpk] = deal(cell(size(csFiles_prm)));
    for iFile=1:numel(csFiles_prm)
        try
            vcFile_prm_ = csFiles_prm{iFile};
            S_score1 = load(strrep(vcFile_prm_, '.prm', '_score.mat'));
            P = loadParam_(vcFile_prm_);
            set0_(P);
            S_ = S_score1.S_score_clu;
            cvrSnr{iFile} = gather_(S_score1.vrSnr_min_gt');
            [cvrFp{iFile}, cvrFn{iFile}, cvrAccuracy{iFile}] = deal(S_.vrFp, S_.vrMiss, S_.vrAccuracy);
            cvnSpk{iFile} = cellfun(@numel, S_.cviSpk_gt_hit) + cellfun(@numel, S_.cviSpk_gt_miss);
            cvnSite{iFile} = S_score1.vnSite_gt;
            disp(csFiles_prm{iFile});
            disp_score_(cvrSnr{iFile}, cvrFp{iFile}, cvrFn{iFile}, cvrAccuracy{iFile}, cvnSite{iFile}, cvnSpk{iFile}, 0);
        catch
            disperr_();
        end
    end

    [vrSnr, vrFp, vrFn, vrAccuracy, vnSite, vnSpk] = multifun_(@(x)cell2mat_(x'), cvrSnr, cvrFp, cvrFn, cvrAccuracy, cvnSite, cvnSpk);
    disp('All files pooled:');
    disp_score_(vrSnr, vrFp, vrFn, vrAccuracy, vnSite, vnSpk, 1);

    % Plot
    vpp_bin = 2; %25 for vpp

    figure;  ax=[];
    ax(1)=subplot(311); hold on;
    plot(vrSnr(:), vrFp(:), 'k.', 'MarkerSize', 2);
    boxplot_(vrFp(:), vrSnr(:), vpp_bin, [3, 40]); ylabel('False Positive');
    xlabel('SNR (Vp/Vrms)'); grid on; set(gca,'YScale','linear');
    ylim_([0, .2]); xlim_([0 40]);
    title_(vcFile_batch);

    ax(2)=subplot(312);  hold on;
    plot(vrSnr(:), vrFn(:), 'k.', 'MarkerSize', 2);
    boxplot_(vrFn(:), vrSnr(:), vpp_bin, [3, 40]); ylabel('False Negative');
    xlabel('SNR (Vp/Vrms)'); grid on; set(gca,'YScale','linear');
    ylim_([0, .2]); xlim_([0 40]);
    set(gcf,'Color','w');

    ax(3)=subplot(313);  hold on;
    plot(vrSnr(:), vnSite(:), 'k.', 'MarkerSize', 2);
    boxplot_(vnSite(:), vrSnr(:), vpp_bin, [3, 40]); ylabel('#sites>thresh');
    xlabel('SNR (Vp/Vrms)'); grid on; set(gca,'YScale','linear');
    linkaxes(ax,'x');
    ylim_([0, .2]); xlim_([0 40]);
    set(gcf,'Color','w');
    ylim_([0, 16]);
end %func
