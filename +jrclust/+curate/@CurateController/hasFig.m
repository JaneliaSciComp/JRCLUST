function flag = hasFig(obj, figKey)
    %HASFIG Return true if we have a figure by key
    flag = ischar(figKey) && isKey(obj.hFigs, figKey);
end