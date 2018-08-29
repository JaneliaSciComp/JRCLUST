%--------------------------------------------------------------------------
function import_gt_silico_(vcFile_mat)
    if matchFileExt_(vcFile_mat, '.mat')
        % catalin's raster function
        S = load(vcFile_mat);
        vnSpk = cellfun(@numel, S.a);
        viClu = int32(cell2mat_(arrayfun(@(n)n*ones(1, vnSpk(n)), 1:numel(vnSpk), 'UniformOutput', 0)));
        viTime = int32(cell2mat_(S.a) * 20); % Convert to sample # (saved in ms unit & sampling rate =20KHZ)
    elseif matchFileExt_(vcFile_mat, '.mda')
        % import Jeremy Magland format
        mr = readmda_(vcFile_mat)';
        viClu = int32(mr(:,3));
        viTime = int32(mr(:,2));
    end
    S_gt = makeStruct(viClu, viTime);
    groundTruthFile = subsFileExt_(vcFile_mat, '_gt.mat');
    write_struct_(groundTruthFile, S_gt);
end % function
