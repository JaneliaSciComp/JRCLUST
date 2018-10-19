%--------------------------------------------------------------------------
% Strip comments from cell string
% 7/24/17 JJJ: Code cleanup
function csLines = strip_comments_(csLines)
    csLines = csLines(cellfun(@(x)~isempty(x), csLines));
    csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0);
    csLines = csLines(cellfun(@(x)x(1)~='%', csLines));

    % remove comments in the middle
    for i=1:numel(csLines)
        vcLine1 = csLines{i};
        iComment = find(vcLine1=='%', 1, 'first');
        if ~isempty(iComment)
            vcLine1 = vcLine1(1:iComment-1);
        end
        vcLine1 = strrep(vcLine1, '...', '');
        if ismember(strsplit(vcLine1), {'for', 'end', 'if'})
            csLines{i} = [strtrim(vcLine1), ', ']; %add blank at the end
        else
            csLines{i} = [strtrim(vcLine1), ' ']; %add blank at the end
        end
    end
    % csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0);
    csLines = csLines(cellfun(@(x)~isempty(x), csLines));
end %func
