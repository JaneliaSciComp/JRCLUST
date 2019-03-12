function closeFilt(obj)
    %CLOSEFILT Close the filtered recording file, clear its data
    if ~obj.filtIsOpen
        return;
    end
    obj.filteredData = [];
    obj.filtIsOpen = 0;
end