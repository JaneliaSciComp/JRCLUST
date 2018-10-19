%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function cellstr2file_(vcFile, csLines)
    % Write a cellstring to a text file

    fid = fopen(vcFile, 'w');
    for i=1:numel(csLines)
        fprintf(fid, '%s\n', csLines{i});
    end
    fclose(fid);
end %func
