%--------------------------------------------------------------------------
function mr = filterq_(vrA, vrB, mr, dimm)
    % quick filter using single instead of filtfilt
    % faster than filtfilt and takes care of the time shift

    if nargin < 4, dimm = 1; end
    if numel(vrA)==1, return; end
    if isempty(vrB), vrB=sum(vrA); end
    %JJJ 2015 09 16
    % mr = filter(vrA, vrB, mr, [], dimm);
    mr = circshift(filter(vrA, vrB, mr, [], dimm), -ceil(numel(vrA)/2), dimm);
end % function
