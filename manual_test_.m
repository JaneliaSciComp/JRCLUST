secondarySelectedCluster%--------------------------------------------------------------------------
function manual_test_(P, csCmd)
    drawnow;
    if nargin<2, csCmd = ''; end
    if isempty(csCmd), csCmd = {'Mouse', 'Menu', 'FigWav', 'FigTime', 'FigClusterCor', 'FigProj', 'Exit'}; end
    if ischar(csCmd), csCmd = {csCmd}; end
    S0 = get(0, 'UserData');
    S_clu = S0.S_clu;
    nClu = S0.S_clu.nClusters;

    for iCmd = 1:numel(csCmd)
        vcCmd1 = csCmd{iCmd};
        fprintf('\tTesting manual-mode %d/%d: %s\n', iCmd, numel(csCmd), vcCmd1);
        switch vcCmd1
            case 'Mouse' % simualte mouse click
            keyPress_fig_(getCachedFig('FigWav'), 'r'); %view whole
            fprintf('\tTesting mouse L/R clicks.\n');
            viClu_test1 = [subsample_vr_(1:nClu, 5), nClu];
            for iClu1=viClu_test1
                fprintf('\t\tprimarySelectedCluster:%d/%d\n', iClu1, numel(viClu_test1));
                update_cursor_([], iClu1, 0);
                keyPressFcn_cell_(getCachedFig('FigWav'), {'c','t','j','i','v','e','f'});
                drawnow;
                viClu_test2 = keep_lim_(iClu1 + [-2:2], [1, nClu]);
                for iClu2=viClu_test2
                    fprintf('\t\t\tsecondarySelectedCluster:%d/%d\n', iClu2, numel(viClu_test2));
                    update_cursor_([], iClu2, 1);
                    keyPressFcn_cell_(getCachedFig('FigWav'), {'c','t','j','i','v','e','f'});
                    drawnow;
                end
            end

            case 'Menu' % run menu items, except for the exit and save (make a black list)
            %             csMenu_skip = {'Show traces', 'Exit'};
            %             hFigWav = figureByTag('FigWav');
            menu_test_(figureByTag('FigWav'), {'Show traces', 'Exit'});
            %             vMenu0 = findobj('Type', 'uimenu', 'Parent', hFigWav);
            %             cvMenu = cell(size(vMenu0));
            %             for iMenu0 = 1:numel(vMenu0)
            %                 cvMenu{iMenu0} = findobj('Type', 'uimenu', 'Parent', vMenu0(iMenu0))';
            %             end
            %             vMenu = [cvMenu{:}];
            %             cCallback_menu = get(vMenu, 'Callback');
            %             csLabel_menu = get(vMenu, 'Label');
            %             fprintf('\tTesting menu items\n');
            %             for iMenu = 1:numel(csLabel_menu)
            %                 vcMenu = csLabel_menu{iMenu};
            %                 if ismember(vcMenu, csMenu_skip), continue; end
            %                 try
            %                     hFunc = cCallback_menu{iMenu};
            % %                     hFunc(hFigWav, []); %call function
            %                     hFunc(vMenu(iMenu), []); %call function
            %                     fprintf('\tMenu ''%s'' success.\n', vcMenu);
            %                 catch
            %                     fprintf(2, '\tMenu ''%s'' failed.\n', vcMenu);
            %                     disperr_();
            %                 end
            %             end

            case 'FigWav' % test all possible keyboard press
            keyPress_fig_(getCachedFig('FigWav'), get_keyPress_('all'));

            case 'FigTime'
            keyPress_fig_(getCachedFig('FigTime'), get_keyPress_('all'));

            case 'FigClusterCor'
            keyPress_fig_(getCachedFig('FigClusterCor'), get_keyPress_('all'));

            case 'FigProj'
            keyPress_fig_(getCachedFig('FigProj'), get_keyPress_('all'));

            case 'Exit'
            %fDebug_ui = 0;  setUserData(fDebug_ui); % disable debug flag
            exit_manual_(getCachedFig('FigWav'));
            %             fDebug_ui = 1;  setUserData(fDebug_ui);

            otherwise
            fprintf(2, 'Unsupported testing mode: %s\n', vcCmd1);
        end %swtich
    end
end %func
