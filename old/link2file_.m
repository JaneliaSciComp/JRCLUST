%--------------------------------------------------------------------------
function csFile = link2file_(csLink)
    csFile = cell(size(csLink));
    for i=1:numel(csLink)
        vcFile1 = csLink{i};
        iBegin = find(vcFile1=='/', 1, 'last'); % strip ?
        if ~isempty(iBegin), vcFile1 = vcFile1(iBegin+1:end); end

        iEnd = find(vcFile1=='?', 1, 'last'); % strip ?
        if ~isempty(iEnd), vcFile1 = vcFile1(1:iEnd-1); end
        csFile{i} = vcFile1;
    end
end %func
