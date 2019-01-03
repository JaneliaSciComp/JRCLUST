function mr = filterq_(vrA, vrB, mr)
    % quick filter using single instead of filtfilt
    % faster than filtfilt and takes care of the time shift

    if numel(vrA) == 1
        return;
    end
    if isempty(vrB)
        vrB = sum(vrA);
    end

    mr = circshift(filter(vrA, vrB, mr, [], 1), -ceil(numel(vrA)/2), 1);
end
