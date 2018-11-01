%--------------------------------------------------------------------------
function keyPressFcn_FigWavCor_(hObject, event)
    S0 = get(0, 'UserData');

    switch lower(event.Key)
        case 'm' % merge
            ui_merge_(S0);

        case 's' % split
            auto_split_(1); % multi

        case 'k' % kilosort template similarity view
            if get_set_(S0.P, 'fImportKsort', 0)
                plot_FigWavCor_(S0, 'simscore');
            end

        case {'d', 'backspace', 'delete'} % delete
            ui_delete_(S0);

        case 'w' % waveform correlation view
            if get_set_(S0.P, 'fImportKsort', 0)
                plot_FigWavCor_(S0, 'wavecor');
            end
    end % switch
end % func
