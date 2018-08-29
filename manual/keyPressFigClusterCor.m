%--------------------------------------------------------------------------
function keyPressFigClusterCor(hObject, event)
    S0 = get(0, 'UserData');

    switch lower(event.Key)
        case 'm' %merge
            ui_merge_(S0);
        case 's' %split
            autoSplit(1); %multi
        case {'d', 'backspace', 'delete'} %delete
            ui_delete_(S0);
        case 'k' % KiloSort template similarity view
            if getOr(S0.P, 'fImportKilosort', 0)
                plotFigClusterCor(S0, 'simscore');
            end
        case 'w' % waveform correlation view
            if getOr(S0.P, 'fImportKilosort', 0)
                plotFigClusterCor(S0, 'wavecor');
            end
    end %switch
end % function
